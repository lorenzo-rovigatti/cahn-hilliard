/*
 * CahnHilliard.cpp
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#include "CahnHilliard.h"

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

	size = N;
	for(int i = 1; i < dims; i++) {
		size *= N;
	}

	rho = RhoMatrix<double>(size, model->N_species());

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
		for(int bin = 0; bin < size; bin++) {
			double modulation = initial_A * std::cos(initial_k * bin);
			for(int i = 0; i < model->N_species(); i++) {
				double average_rho = densities[i];
				rho(bin, i) = average_rho * (1.0 + modulation * (1.0 + 2.0 * (drand48() - 0.5) * 1e-2));
			}
		}
	}

	dx *= _user_to_internal; // proportional to m
	M /= _user_to_internal; // proportional to m^-1
	k_laplacian *= SQR(SQR(_user_to_internal)) * _user_to_internal; // proportional to m^5
	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho(idx, species) /= CUB(_user_to_internal); // proportional to m^-3
		}
	}

	if(_reciprocal) {
		rho_hat.resize(size / 2 + 1);
		sqr_wave_vectors.resize(size / 2 + 1);
		dealiaser.resize(size / 2 + 1);

		// rho_plan = fftw_plan_dft_r2c_1d(size, rho.data(), reinterpret_cast<fftw_complex *>(rho_hat.data()), FFTW_ESTIMATE);
		// // c2r transforms overwrite the input if FFTW_PRESERVE_INPUT is not specified
		// rho_inverse_plan = fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex *>(rho_hat.data()), rho.data(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

		double nyquist_mode = size * M_PI / (N * dx) * 2.0 / 3.0;
		for(unsigned int i = 0; i < sqr_wave_vectors.size(); i++) {
			double k = 2.0 * M_PI * i / (N * dx);
			dealiaser[i] = (k < nyquist_mode) ? 1.0 : 0.0;
			sqr_wave_vectors[i] = SQR(k);
		}

		fftw_execute(rho_plan);
	}

	_init_CUDA(config);
}

template<int dims>
CahnHilliard<dims>::~CahnHilliard() {
#ifndef NOCUDA
	if(_d_rho != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho));
	}
	if(_d_rho_der != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho_der));
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
		model->der_bulk_free_energy(_d_rho, _d_rho_der, rho.bins());
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
				// if(species == 1 && idx == 0) printf("0 %lf\n", model->der_bulk_free_energy(species, rho[idx]));
				rho_der(idx, species) = model->der_bulk_free_energy(species, rho.rho_species(idx)) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
			}
		}

		// and then we integrate them
		for(unsigned int idx = 0; idx < rho.bins(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				// if(species == 1 && idx == 0) printf("0 %e %lf\n", cell_laplacian(rho_der, species, idx), rho_der[idx][species]);
				rho(idx, species) += M * cell_laplacian(rho_der, species, idx) * dt;
			}
		}
	}
}

template<int dims>
void CahnHilliard<dims>::_evolve_reciprocal() {
	
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	if(!_output_ready) _GPU_CPU();

	double V_bin = 1;
	for(int d = 0; d < dims; d++) {
		V_bin *= dx;
	}

	double mass = 0.;
	for(unsigned int i = 0; i < rho.bins(); i++) {
		mass += rho.rho_tot(i);
	}

	return mass * V_bin;
}

template<int dims>
double CahnHilliard<dims>::total_free_energy() {
	if(!_output_ready) _GPU_CPU();

	double V_bin = 1;
	for(int d = 0; d < dims; d++) {
		V_bin *= dx;
	}

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

	return fe * V_bin;
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

	for(int idx = 0; idx < size; idx++) {
		output << _density_to_user(rho(idx, species)) << " " << std::endl;
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, std::ofstream &output) {
	if(!_output_ready) _GPU_CPU();

	for(int idx = 0; idx < size; idx++) {
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

	for(int idx = 0; idx < size; idx++) {
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
	CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_der, d_der_vec_size)); // float instead of double

	init_symbols(N, size, model->N_species());

	_CPU_GPU();
#endif
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
