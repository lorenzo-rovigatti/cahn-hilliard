/*
 * CahnHilliard.cpp
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#include "CahnHilliard.h"

#include "integrators/BailoFiniteVolume.h"
#include "integrators/EulerCPU.h"
#include "integrators/EulerMobilityCPU.h"
#include "integrators/PseudospectralCPU.h"

#include "utils/utility_functions.h"

#ifndef NOCUDA
#include "CUDA/integrators/EulerCUDA.h"
#include "CUDA/integrators/EulerMobilityCUDA.h"
#include "CUDA/integrators/PseudospectralCUDA.h"
#endif

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
	dx = _config_optional_value<double>(config, "dx", 1.0);
	_internal_to_user = _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
	_user_to_internal = 1.0 / _internal_to_user;
	std::string user_integrator = _config_optional_value<std::string>(config, "integrator", "euler");

	bool use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);

	info("Running a simulation with N = {}, dt = {}, dx = {}, scaling factor = {}", N, dt, dx, _internal_to_user);

	double log2N = std::log2(N);
	if(ceil(log2N) != floor(log2N)) {
		critical("N should be a power of 2");
	}

	N_minus_one = N - 1;
	bits = (int) log2N;

	grid_size = N;
	_grid_size_str = fmt::format("{}", N);
	for(int i = 1; i < dims; i++) {
		_grid_size_str += fmt::format("x{}", N);
		grid_size *= N;
	}

	RhoMatrix<double> rho(grid_size, model->N_species());

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
		for(int bin = 0; bin < grid_size; bin++) {
			double modulation = initial_A * std::cos(initial_k * bin);
			for(int i = 0; i < model->N_species(); i++) {
				double random_factor = (initial_N_peaks == 0) ? (drand48() - 0.5) : 1.0 + 0.02 * (drand48() - 0.5);
				double average_rho = densities[i];
				rho(bin, i) = average_rho * (1.0 + 2.0 * modulation * random_factor);
			}
		}
	}

	dx *= _user_to_internal; // proportional to m
	k_laplacian *= std::pow(_user_to_internal, 5); // proportional to m^5
	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho(idx, species) /= CUB(_user_to_internal); // proportional to m^-3
		}
	}

	V_bin = CUB(dx);

	if(user_integrator == "euler") {
		std::string mobility = _config_optional_value<std::string>(config, "mobility.type", "constant");

		if(use_CUDA) {
#ifndef NOCUDA
			if(mobility != "constant") {
				integrator = new EulerMobilityCUDA<dims>(m, config);
			}
			else {
				integrator = new EulerCUDA<dims>(m, config);
			}
#endif
		}
		else {
			if(mobility == "constant") {
				integrator = new EulerCPU<dims>(m, config);
			}
			else {
				integrator = new EulerMobilityCPU<dims>(m, config);
			}
		}
	}
	else if(user_integrator == "pseudospectral") {
		if(use_CUDA) {
#ifndef NOCUDA
			integrator = new PseudospectralCUDA<dims>(m, config);
#endif
		}
		else {
			integrator = new PseudospectralCPU<dims>(m, config);
		}
	}
	else if(user_integrator == "bailo") {
		if(use_CUDA) {
#ifndef NOCUDA
			critical("Unsupported CUDA integrator {}", user_integrator);
#endif
		}
		else {
			integrator = new BailoFiniteVolume<dims>(m, config);
		}
	}
	else {
		critical("Unsupported integrator {}", user_integrator);
	}
	
	integrator->set_initial_rho(rho);
}

template<int dims>
CahnHilliard<dims>::~CahnHilliard() {
	if(integrator != nullptr) {
		delete integrator;
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

template<int dims>
void CahnHilliard<dims>::evolve() {
	integrator->evolve();
}

template<int dims>
double CahnHilliard<dims>::average_mass() {
	double mass = 0.;
	for(unsigned int i = 0; i < integrator->rho().bins(); i++) {
		mass += integrator->rho().rho_tot(i);
		if(safe_isnan(mass)) {
			critical("Encountered a nan while computing the total mass (bin {})", i);
		}
	}

	return mass * V_bin / integrator->rho().bins();
}

template<int dims>
double CahnHilliard<dims>::average_free_energy() {
	double fe = 0.;
	for(unsigned int i = 0; i < integrator->rho().bins(); i++) {
		double interfacial_contrib = 0.;
		for(int species = 0; species < model->N_species(); species++) {
			auto rho_grad = gradient(integrator->rho(), species, i);
			for(int d = 0; d < dims; d++) {
				interfacial_contrib += k_laplacian * rho_grad[d] * rho_grad[d];
			}
		}
		fe += model->bulk_free_energy(integrator->rho().rho_species(i)) + interfacial_contrib;
	}

	return fe * V_bin / integrator->rho().bins();
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, const std::string &filename, long long int t) {
	std::ofstream output(filename);
	print_species_density(species, output, t);
	output.close();
}

template<>
void CahnHilliard<1>::print_species_density(int species, std::ofstream &output, long long int t) {
	output << fmt::format("# step = {}, t = {:.5}, size = {}", t, t * dt, _grid_size_str) << std::endl;
	for(int idx = 0; idx < grid_size; idx++) {
		output << _density_to_user(integrator->rho()(idx, species)) << " " << std::endl;
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, std::ofstream &output, long long int t) {
	output << fmt::format("# step = {}, t = {:.5}, size = {}", t, t * dt, _grid_size_str) << std::endl;
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
		output << _density_to_user(integrator->rho()(idx, species)) << " ";
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_total_density(const std::string &filename, long long int t) {
	std::ofstream output(filename.c_str());

	output << fmt::format("# step = {}, t = {:.5}, size = {}", t, t * dt, _grid_size_str) << std::endl;
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
		output << _density_to_user(integrator->rho().rho_tot(idx)) << std::endl;
	}

	output.close();
}


template<int dims>
double CahnHilliard<dims>::_density_to_user(double v) {
	return v / (_internal_to_user * _internal_to_user * _internal_to_user);
}

template class CahnHilliard<1>;
template class CahnHilliard<2>;

} /* namespace ch */
