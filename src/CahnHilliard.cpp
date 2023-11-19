/*
 * CahnHilliard.cpp
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#include "CahnHilliard.h"

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

	info("Running a simulation with N = {}, dt = {}, H = {}, M = {}", N, dt, dx, M);

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

	rho.resize(size, std::vector<double>(model->N_species(), 0.));

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
					load_from >> rho[idx][s];
				}
				break;
			case 2:
				int coords[2];
				for(coords[1] = 0; coords[1] < N; coords[1]++) {
					for(coords[0] = 0; coords[0] < N; coords[0]++) {
						int idx = cell_idx(coords);
						load_from >> rho[idx][s];
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
		for(int bin = 0; bin < N; bin++) {
			auto &species_rho = rho[bin];
			double modulation = initial_A * std::cos(initial_k * bin);
			for(int i = 0; i < species_rho.size(); i++) {
				double average_rho = densities[i];
				species_rho[i] = average_rho * (1.0 + modulation * (1.0 + 2.0 * (drand48() - 0.5) * 1e-2));
			}
		}
	}

	_init_CUDA(config);
}

template<int dims>
CahnHilliard<dims>::~CahnHilliard() {
	if(_d_rho != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho));
	}
	if(_d_rho_der != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho_der));
	}
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
double CahnHilliard<1>::cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx) {
	int idx_m = (idx - 1 + N) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (field[idx_m][species] + field[idx_p][species] - 2.0 * field[idx][species]) / SQR(dx);
}

template<>
double CahnHilliard<2>::cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx) {
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
			field[cell_idx(coords_xmy)][species] +
			field[cell_idx(coords_xpy)][species] +
			field[cell_idx(coords_xym)][species] +
			field[cell_idx(coords_xyp)][species] -
			4 * field[idx][species])
			/ SQR(dx);
}

template<int dims>
void CahnHilliard<dims>::evolve() {
	if(_use_CUDA) {
		model->der_bulk_free_energy(_d_rho, _d_rho_der, rho.size());
		_output_ready = false;
	}
	else {
		static std::vector<std::vector<double>> rho_der(rho.size(), std::vector<double>(model->N_species()));
		// we first evaluate the time derivative for all the fields
		for(unsigned int idx = 0; idx < rho.size(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				rho_der[idx][species] = model->der_bulk_free_energy(species, rho[idx]) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
			}
		}

		// and then we integrate them
		for(unsigned int idx = 0; idx < rho.size(); idx++) {
			for(int species = 0; species < model->N_species(); species++) {
				rho[idx][species] += M * cell_laplacian(rho_der, species, idx) * dt;
			}
		}
	}
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	if(!_output_ready) _GPU_CPU();

	double mass = 0.;
	for(unsigned int i = 0; i < rho.size(); i++) {
		mass += std::accumulate(rho[i].begin(), rho[i].end(), 0.);
	}

	return mass;
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
		output << rho[idx][species] << " " << std::endl;
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
		output << rho[idx][species] << " ";
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
		output << std::accumulate(rho[idx].begin(), rho[idx].end(), 0.) << std::endl;
	}

	output.close();
}

template<int dims>
void CahnHilliard<dims>::_init_CUDA(toml::table &config) {
	_use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);
	if(!_use_CUDA) return;

	_d_vec_size = rho.size() * model->N_species() * sizeof(number);

	info("Initialising CUDA arrays of size {} bytes", _d_vec_size);

	_h_rho.resize(rho.size() * model->N_species());
	CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho, _d_vec_size));
	CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_der, _d_vec_size));

	_CPU_GPU();
}

template<int dims>
void CahnHilliard<dims>::_CPU_GPU() {
	if(!_use_CUDA) return;

	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			_h_rho[rho.size() * species + idx] = rho[idx][species];
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_rho, _h_rho.data(), _d_vec_size, cudaMemcpyHostToDevice));
}

template<int dims>
void CahnHilliard<dims>::_GPU_CPU() {
	if(!_use_CUDA) return;

	CUDA_SAFE_CALL(cudaMemcpy(_h_rho.data(), _d_rho, _d_vec_size, cudaMemcpyDeviceToHost));

	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho[idx][species] = _h_rho[rho.size() * species + idx];
		}
	}

	_output_ready = true;
}

template class CahnHilliard<1>;
template class CahnHilliard<2>;

} /* namespace ch */