/*
 * CahnHilliard.cpp
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#include "CahnHilliard.h"

#include "utils/utility_functions.h"
#include "CUDA/CahnHilliard.cuh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace ch {

template<int dims>
CahnHilliard<dims>::CahnHilliard(FreeEnergyModel *m, toml::table &config) :
				model(m) {

	N = _config_value<int>(config, "N");
	k_laplacian = _config_optional_value<double>(config, "k", 1.0);
	dt = _config_value<double>(config, "dt");
	M = _config_optional_value<double>(config, "M", 1.0);
	dx = _config_optional_value<double>(config, "dx", 1.0);
	_internal_to_user = _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
	_user_to_internal = 1.0 / _internal_to_user;

	_reciprocal = _config_optional_value<bool>(config, "use_reciprocal_space", false);

	info("Running a simulation with N = {}, dt = {}, dx = {}, M = {}, scaling factor = {}", N, dt, dx, M, _internal_to_user);

	double log2N = std::log2(N);
	if(ceil(log2N) != floor(log2N)) {
		critical("N should be a power of 2");
	}

	N_minus_one = N - 1;
	bits = (int) log2N;

	grid_size = N;
	for(int i = 1; i < dims; i++) {
		grid_size *= N;
	}

	rho = RhoMatrix<double>(grid_size, model->N_species());

	if(!config["initial_density"] && !config["load_from"]) {
		critical("Either 'initial_density' or 'load_from' should be specified");
	}

	if(config["load_from"]) {
		std::string filename = _config_value<std::string>(config, "load_from");
		std::ifstream load_from(filename.c_str());

		if(load_from.good()) {
			info("Using filename '{}' for the initial conditions", filename);
		}
		else {
			critical("The initial conditions file '{}' is not readable", filename);
		}

		for(int s = 0; s < model->N_species(); s++) {
			switch(dims) {
			case 1:
				for(int idx = 0; idx < N; idx++) {
					load_from >> rho(idx, s);
				}
				break;
			case 2:
				int coords[2];
				for(coords[1] = 0; coords[1] < N; coords[1]++) {
					for(coords[0] = 0; coords[0] < N; coords[0]++) {
						int idx = cell_idx(coords);
						load_from >> rho(idx, s);
					}
				}

				break;
			default:
				critical("Unsupported number of dimensions {}", dims);
			}
		}

		load_from.close();
	}
	else { // initial-density
		std::vector<double> densities = _config_array_values<double>(config, "initial_density", model->N_species());
		double initial_A = _config_optional_value<double>(config, "initial_A", 1e-2);
		int initial_N_peaks = _config_optional_value<int>(config, "initial_N_peaks", 0);
		double initial_k = 2 * M_PI * initial_N_peaks / (double) N; // wave vector of the modulation
		int seed = _config_optional_value<int>(config, "seed", time(NULL));
		srand48(seed);
		for(int bin = 0; bin < grid_size; bin++) {
			double modulation = initial_A * std::cos(initial_k * bin);
			for(int i = 0; i < model->N_species(); i++) {
				double average_rho = densities[i];
				rho(bin, i) = average_rho * (1.0 + modulation * (1.0 + 2.0 * (drand48() - 0.5) * 1e-2));
			}
		}
	}

	dx *= _user_to_internal; // proportional to m
	M /= _user_to_internal; // proportional to m^-1
	k_laplacian *= std::pow(_user_to_internal, 5); // proportional to m^5
	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho(idx, species) /= CUB(_user_to_internal); // proportional to m^-3
		}
	}

	V_bin = CUB(dx);

	if(_reciprocal) {
		_reciprocal_n.fill(N);
		hat_grid_size = _reciprocal_n[dims - 1] / 2 + 1; 
		for(int i = 0; i < dims - 1; i++) {
			hat_grid_size *= _reciprocal_n[i];
		}
		int hat_vector_size = hat_grid_size * model->N_species();

		info("Size of the reciprocal vectors: {}", hat_vector_size);

		rho_hat.resize(hat_vector_size);
		rho_hat_copy.resize(hat_vector_size);
		sqr_wave_vectors.resize(hat_vector_size);
		dealiaser.resize(hat_vector_size);
		f_der = RhoMatrix<double>(rho.bins(), model->N_species());
		f_der_hat.resize(hat_vector_size);
		
		int idist = grid_size;
		int odist = hat_grid_size; // the distance between the first elements of arrays referring to nearby species in the reciprocal space
		f_der_plan = fftw_plan_many_dft_r2c(dims, _reciprocal_n.data(), model->N_species(), f_der.data(), NULL, 1, idist, reinterpret_cast<fftw_complex *>(f_der_hat.data()), NULL, 1, odist, FFTW_ESTIMATE);

		// c2r transforms overwrite the input array
		rho_inverse_plan = fftw_plan_many_dft_c2r(dims, _reciprocal_n.data(), model->N_species(), reinterpret_cast<fftw_complex *>(rho_hat_copy.data()), NULL, 1, odist, rho.data(), NULL, 1, idist, FFTW_ESTIMATE);

		double nyquist_mode = N * M_PI / (N * dx) * 2.0 / 3.0;
		if(dims == 1) {
			int k_idx = 0;
			for(int species = 0; species < model->N_species(); species++) {
				for(unsigned int i = 0; i < hat_grid_size; i++) {
					double k = 2.0 * M_PI * i / (N * dx);
					dealiaser[k_idx] = (k < nyquist_mode) ? 1.0 : 0.0;
					sqr_wave_vectors[k_idx] = SQR(k);
					k_idx++;
				}
			}
		}
		else if(dims == 2) {
			int k_idx = 0;
			for(int species = 0; species < model->N_species(); species++) {
				for(int kx_idx = 0; kx_idx < _reciprocal_n[0]; kx_idx++) {
					int kx = (kx_idx < _reciprocal_n[0] / 2) ? kx_idx : _reciprocal_n[0] - kx_idx;
					for(int ky = 0; ky < (_reciprocal_n[1] / 2 + 1); ky++) {
						double k = 2.0 * M_PI * std::sqrt(SQR(kx) + SQR(ky)) / (N * dx);
						dealiaser[k_idx] = (k < nyquist_mode) ? 1.0 : 0.0;
						sqr_wave_vectors[k_idx] = SQR(k);
						k_idx++;
					}
				}
			}
		}
		else {
			critical("Unsupported number of dimensions {}", dims);
		}

		fftw_plan rho_plan = fftw_plan_many_dft_r2c(dims, _reciprocal_n.data(), model->N_species(), rho.data(), NULL, 1, idist, reinterpret_cast<fftw_complex *>(rho_hat.data()), NULL, 1, odist, FFTW_ESTIMATE);
		fftw_execute(rho_plan);
		rho_hat_copy = rho_hat;
		fftw_destroy_plan(rho_plan);
	}

	_init_CUDA(config);
}

template<int dims>
CahnHilliard<dims>::~CahnHilliard() {
	if(_reciprocal) {
		fftw_destroy_plan(rho_inverse_plan);
		fftw_destroy_plan(f_der_plan);
		fftw_cleanup();
	}

#ifndef NOCUDA
	if(_d_rho != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho));
	}
	if(_d_rho_der != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho_der));
	}
	if(_d_rho_hat != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho_hat));
		CUDA_SAFE_CALL(cudaFree(_d_rho_hat_copy));
		CUFFT_CALL(cufftDestroy(_d_rho_inverse_plan));
	}
	if(_d_f_der_hat != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_f_der_hat));
		CUFFT_CALL(cufftDestroy(_d_f_der_plan));
	}
	if(_d_sqr_wave_vectors != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_sqr_wave_vectors));
	}
	if(_d_dealiaser != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_dealiaser));
	}
#endif
}

template<int dims>
void CahnHilliard<dims>::fill_coords(int coords[dims], int idx) {
	for(int d = 0; d < dims; d++) {
		coords[d] = idx & N_minus_one;
		idx >>= bits; // divide by N
	}
}

template<int dims>
int CahnHilliard<dims>::cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= bits; // multiply by N
	}
	return idx;
}

template<>
std::array<double, 1> CahnHilliard<1>::gradient(RhoMatrix<double> &field, int species, int idx) {
	int idx_p = (idx + 1) & N_minus_one;

	return {(field(idx_p, species) - field(idx, species)) / dx};
}

template<>
std::array<double, 2> CahnHilliard<2>::gradient(RhoMatrix<double> &field, int species, int idx) {
	int coords_xy[2];
	fill_coords(coords_xy, idx);

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & N_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & N_minus_one
	};

	return {
		(field(cell_idx(coords_xpy), species) - field(idx, species)) / dx,
		(field(cell_idx(coords_xyp), species) - field(idx, species)) / dx
	};
}

template<>
double CahnHilliard<1>::cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + N) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(dx);
}

template<>
double CahnHilliard<2>::cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
	int coords_xy[2];
	fill_coords(coords_xy, idx);

	int coords_xmy[2] = {
			(coords_xy[0] - 1 + N) & N_minus_one,
			coords_xy[1]
	};

	int coords_xym[2] = {
			coords_xy[0],
			(coords_xy[1] - 1 + N) & N_minus_one
	};

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & N_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & N_minus_one
	};

	return (
			field(cell_idx(coords_xmy), species) +
			field(cell_idx(coords_xpy), species) +
			field(cell_idx(coords_xym), species) +
			field(cell_idx(coords_xyp), species) -
			4 * field(idx, species))
			/ SQR(dx);
}

template<int dims>
void CahnHilliard<dims>::evolve() {
	if(_reciprocal) {
		_evolve_reciprocal();
	}
	else {
		_evolve_direct();
	}
}

template<int dims>
void CahnHilliard<dims>::_evolve_direct() {
	if(_use_CUDA) {
#ifndef NOCUDA
		model->der_bulk_free_energy(_d_rho, _d_rho_der, model->N_species() * rho.bins());
		add_surface_term<dims>(_d_rho, _d_rho_der, dx, k_laplacian);
		integrate<dims>(_d_rho, _d_rho_der, dx, dt, M);

		_output_ready = false;
#endif
	}
	else {
		static RhoMatrix<double> rho_der(rho.bins(), model->N_species());
		// we first evaluate the time derivative for all the fields
		for(unsigned int idx = 0; idx < rho.bins(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				rho_der(idx, species) = model->der_bulk_free_energy(species, rho.rho_species(idx)) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
			}
		}

		// and then we integrate them
		for(unsigned int idx = 0; idx < rho.bins(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				rho(idx, species) += M * cell_laplacian(rho_der, species, idx) * dt;
			}
		}
	}
}

template<int dims>
void CahnHilliard<dims>::_evolve_reciprocal() {
	if(_use_CUDA) {
#ifndef NOCUDA
		model->der_bulk_free_energy(_d_rho, _d_rho_der, rho.bins());

		CUFFT_CALL(cufftExecR2C(_d_f_der_plan, _d_rho_der, _d_f_der_hat));

		integrate_fft(_d_rho_hat, _d_rho_hat_copy, _d_f_der_hat, _d_sqr_wave_vectors, _d_dealiaser, dt, M, k_laplacian, hat_vector_size);

#ifdef CUDA_FIELD_FLOAT
		CUFFT_CALL(cufftExecC2R(_d_rho_inverse_plan, _d_rho_hat_copy, _d_rho));
#else
		CUFFT_CALL(cufftExecZ2D(_d_rho_inverse_plan, _d_rho_hat_copy, _d_rho));
#endif

		_output_ready = false;
#endif
	}
	else {
		for(unsigned int idx = 0; idx < rho.bins(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				f_der(idx, species) = model->der_bulk_free_energy(species, rho.rho_species(idx));
			}
		}

		fftw_execute(f_der_plan); // transform f_der into f_der_hat

		for(unsigned int k_idx = 0; k_idx < rho_hat.size(); k_idx++) {
			// f_der_hat[k_idx] *= dealiaser[k_idx];
			rho_hat[k_idx] = rho_hat_copy[k_idx] = (rho_hat[k_idx] - dt * M * sqr_wave_vectors[k_idx] * f_der_hat[k_idx]) / (1.0 + dt * M * 2.0 * k_laplacian * SQR(sqr_wave_vectors[k_idx]));
			rho_hat_copy[k_idx] /= grid_size;
		}

		fftw_execute(rho_inverse_plan);
	}
}

template<int dims>
double CahnHilliard<dims>::average_mass() {
	if(!_output_ready) _GPU_CPU();

	double mass = 0.;
	for(unsigned int i = 0; i < rho.bins(); i++) {
		mass += rho.rho_tot(i);
		if(safe_isnan(mass)) {
			critical("Encountered a nan while computing the total mass (bin {})", i);
		}
	}

	return mass * V_bin / rho.bins();
}

template<int dims>
double CahnHilliard<dims>::average_free_energy() {
	if(!_output_ready) _GPU_CPU();

	double fe = 0.;
	for(unsigned int i = 0; i < rho.bins(); i++) {
		double interfacial_contrib = 0.;
		for(int species = 0; species < model->N_species(); species++) {
			auto rho_grad = gradient(rho, species, i);
			for(int d = 0; d < dims; d++) {
				interfacial_contrib += k_laplacian * rho_grad[d] * rho_grad[d];
			}
		}
		fe += model->bulk_free_energy(rho.rho_species(i)) + interfacial_contrib;
	}

	return fe * V_bin / rho.bins();
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, const std::string &filename) {
	if(!_output_ready) _GPU_CPU();

	std::ofstream output(filename);
	print_species_density(species, output);
	output.close();
}

template<>
void CahnHilliard<1>::print_species_density(int species, std::ofstream &output) {
	if(!_output_ready) _GPU_CPU();

	for(int idx = 0; idx < grid_size; idx++) {
		output << _density_to_user(rho(idx, species)) << " " << std::endl;
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, std::ofstream &output) {
	if(!_output_ready) _GPU_CPU();

	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0) {
			int modulo = N;
			for(int d = 1; d < dims; d++) {
				if(idx % modulo == 0) {
					output << std::endl;
				}
				modulo <<= bits;
			}
		}
		output << _density_to_user(rho(idx, species)) << " ";
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_total_density(const std::string &filename) {
	if(!_output_ready) _GPU_CPU();

	std::ofstream output(filename.c_str());

	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0) {
			int modulo = N;
			for(int d = 1; d < dims; d++) {
				if(idx % modulo == 0) {
					output << std::endl;
				}
				modulo <<= bits;
			}
		}
		output << _density_to_user(rho.rho_tot(idx)) << std::endl;
	}

	output.close();
}


template<int dims>
double CahnHilliard<dims>::_density_to_user(double v) {
	return v / (_internal_to_user * _internal_to_user * _internal_to_user);
}

template<int dims>
void CahnHilliard<dims>::_init_CUDA(toml::table &config) {
#ifndef NOCUDA
	_use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);
	if(!_use_CUDA) return;

	_d_vec_size = rho.bins() * model->N_species() * sizeof(field_type);
	int d_der_vec_size = rho.bins() * model->N_species() * sizeof(float);

	info("Initialising CUDA arrays of size {} ({} bytes)", rho.bins() * model->N_species(), _d_vec_size);

	_h_rho = RhoMatrix<field_type>(rho.bins(), model->N_species());
	CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho, _d_vec_size));
	CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_der, d_der_vec_size)); // always float

	init_symbols(N, grid_size, model->N_species());

	if(_reciprocal) {
		CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_hat, sizeof(cufftFieldComplex) * rho_hat.size()));
		CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_hat_copy, sizeof(cufftFieldComplex) * rho_hat.size()));
		CUDA_SAFE_CALL(cudaMalloc((void **) &_d_f_der_hat, sizeof(cufftComplex) * f_der_hat.size()));
		CUDA_SAFE_CALL(cudaMalloc((void **) &_d_sqr_wave_vectors, sizeof(float) * sqr_wave_vectors.size()));
		CUDA_SAFE_CALL(cudaMalloc((void **) &_d_dealiaser, sizeof(float) * dealiaser.size()));

		// sqr_wave_vectors e deliaser are std::vector<double>, so we first have to convert them to std::vector<float>
		// and then we can copy their content to the GPU
		std::vector<float> f_sqr_wave_vectors(sqr_wave_vectors.begin(), sqr_wave_vectors.end());
		CUDA_SAFE_CALL(cudaMemcpy(_d_sqr_wave_vectors, f_sqr_wave_vectors.data(), sizeof(float) * f_sqr_wave_vectors.size(), cudaMemcpyHostToDevice));
		std::vector<float> f_dealiaser(dealiaser.begin(), dealiaser.end());
		CUDA_SAFE_CALL(cudaMemcpy(_d_dealiaser, f_dealiaser.data(), sizeof(float) * f_dealiaser.size(), cudaMemcpyHostToDevice));

		CUFFT_CALL(cufftCreate(&_d_rho_inverse_plan));
		CUFFT_CALL(cufftCreate(&_d_f_der_plan));

#ifdef CUDA_FIELD_FLOAT
		CUFFT_CALL(cufftPlanMany(&_d_rho_inverse_plan, dims, _reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_C2R, model->N_species()));
#else
		CUFFT_CALL(cufftPlanMany(&_d_rho_inverse_plan, dims, _reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_Z2D, model->N_species()));
#endif 
		CUFFT_CALL(cufftPlanMany(&_d_f_der_plan, dims, _reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_R2C, model->N_species()));
		
		// copy the initial conditions
		CUDA_SAFE_CALL(cudaMemcpy(_d_rho_hat, rho_hat.data(), sizeof(cufftFieldComplex) * rho_hat.size(), cudaMemcpyHostToDevice));
	}

	_CPU_GPU();
#endif /* NOCUDA */
}

template<int dims>
void CahnHilliard<dims>::_CPU_GPU() {
#ifndef NOCUDA
	if(!_use_CUDA) return;

	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			_h_rho(idx, species) = rho(idx, species);
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_rho, _h_rho.data(), _d_vec_size, cudaMemcpyHostToDevice));
#endif
}

template<int dims>
void CahnHilliard<dims>::_GPU_CPU() {
#ifndef NOCUDA
	if(!_use_CUDA) return;

	CUDA_SAFE_CALL(cudaMemcpy(_h_rho.data(), _d_rho, _d_vec_size, cudaMemcpyDeviceToHost));

	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho(idx, species) = _h_rho(idx, species);
		}
	}

	_output_ready = true;
#endif
}

template class CahnHilliard<1>;
template class CahnHilliard<2>;

} /* namespace ch */
