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
#include "integrators/PseudospectralMobilityCPU.h"

#include "utils/utility_functions.h"
#include "utils/strings.h"

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
CahnHilliard<dims>::CahnHilliard(SimulationState<dims> &sim_state, FreeEnergyModel *m, toml::table &config) :
				_sim_state(sim_state),
				model(m) {

	N = _config_value<int>(config, "N");
	k_laplacian = _config_optional_value<double>(config, "k", 1.0);
	dt = _config_value<double>(config, "dt");
	dx = _config_optional_value<double>(config, "dx", 1.0);
	_internal_to_user = _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
	_user_to_internal = 1.0 / _internal_to_user;
	std::string user_integrator = _config_optional_value<std::string>(config, "integrator", "euler");

	sim_state.use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);

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

	_sim_state.rho = MultiField<double>(grid_size, model->N_species());
	_sim_state.mobility = MultiField<double>(grid_size, model->N_species());

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

		std::string line;
		int lines = 0;
		for(int s = 0; s < model->N_species(); s++) {
			int OK_lines = 0;
			while(OK_lines < N && load_from.good()) {
				std::getline(load_from, line);
				utils::trim(line);
				auto spl = utils::split(line);

				if(line.size() > 0 && line[0] != '#') {
					if(dims == 1) {
						if(spl.size() != 1) {
							critical("Line n. {} in the initial configuration file contains {} fields, should be 1", lines, spl.size());
						}
						_sim_state.rho(OK_lines, s) = std::stod(line);
					}
					else if(dims == 2) {
						if(spl.size() != N) {
							critical("Line n. {} in the initial configuration file contains only {} fields, should be {}", lines, spl.size(), N);
						}
						int coords[2] = {0, OK_lines};
						for(coords[0] = 0; coords[0] < N; coords[0]++) {
							int idx = cell_idx(coords);
							_sim_state.rho(idx, s) = std::stod(spl[coords[0]]);
						}
					}
					else {
						static_assert(dims == 1 || dims == 2, "Unsupported dims for the initial configuration");
					}
					OK_lines++;
				}
				lines++;
			}

			if(OK_lines != N) {
				critical("The initial configuration file contains only {} valid lines for species {}, should be {}", OK_lines, s, N);
			}
		}

		info("Initial configuration parsed: found {} lines, of which {} with data", lines, model->N_species() * N);

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
				if(average_rho != 0.0) {
					_sim_state.rho(bin, i) = average_rho * (1.0 + 2.0 * modulation * random_factor);
				}
				else {
					_sim_state.rho(bin, i) = 2.0 * modulation * random_factor;
				}
			}
		}
	}

	dx *= _user_to_internal; // proportional to m
	k_laplacian *= std::pow(_user_to_internal, 5); // proportional to m^5
	for(unsigned int idx = 0; idx < _sim_state.rho.bins(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			_sim_state.rho(idx, species) /= CUB(_user_to_internal); // proportional to m^-3
		}
	}

	V_bin = CUB(dx);

	// we need to build the mobility object before building the integrator, since CUDA integrators
	// need to allocate device memory for mobility
	mobility = std::unique_ptr<IMobility<dims>>(build_mobility(config, _sim_state));

	if(user_integrator == "euler") {
		if(sim_state.use_CUDA) {
#ifndef NOCUDA
			integrator = new EulerCUDA<dims>(_sim_state, m, config);
#endif
		}
		else {
			integrator = new EulerCPU<dims>(_sim_state, m, config);
		}
	}
	else if(user_integrator == "euler_mobility") {
		if(sim_state.use_CUDA) {
#ifndef NOCUDA
			integrator = new EulerMobilityCUDA<dims>(_sim_state, m, config);
#endif
		}
		else {
			integrator = new EulerMobilityCPU<dims>(_sim_state, m, config);
		}
	}
	else if(user_integrator == "pseudospectral") {
		if(sim_state.use_CUDA) {
#ifndef NOCUDA
			integrator = new PseudospectralCUDA<dims>(_sim_state, m, config);
#endif
		}
		else {
			integrator = new PseudospectralCPU<dims>(_sim_state, m, config);
		}
	}
	else if(user_integrator == "bailo") {
		if(sim_state.use_CUDA) {
#ifndef NOCUDA
			critical("Unsupported CUDA integrator {}", user_integrator);
#endif
		}
		else {
			integrator = new BailoFiniteVolume<dims>(_sim_state, m, config);
		}
	}
	else if(user_integrator == "pseudospectral_mobility") {
		if(sim_state.use_CUDA) {
#ifndef NOCUDA
			critical("Unsupported CUDA integrator {}", user_integrator);
#endif
		}
		else {
			integrator = new PseudospectralMobilityCPU<dims>(_sim_state, m, config);
		}
	}
	else {
		critical("Unsupported integrator {}", user_integrator);
	}
	integrator->validate();

	// the "free_energy" mobility model requires some of the model's internals to be initialised
	// and since mobility is computed before the first evolve() call, we need to set it up now
	MultiField<double> dummy_rho_der(grid_size, model->N_species());
	model->der_bulk_free_energy(_sim_state.rho, dummy_rho_der); // dummy call to initialise the model
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
std::array<double, 1> CahnHilliard<1>::gradient(MultiField<double> &field, int species, int idx) {
	int idx_p = (idx + 1) & N_minus_one;

	return {(field(idx_p, species) - field(idx, species)) / dx};
}

template<>
std::array<double, 2> CahnHilliard<2>::gradient(MultiField<double> &field, int species, int idx) {
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
	mobility->update_mobility();
	integrator->evolve();
}

template<int dims>
double CahnHilliard<dims>::average_mass() {
	integrator->sync();

	double mass = 0.;
	for(unsigned int i = 0; i < _sim_state.rho.bins(); i++) {
		mass += _sim_state.rho.field_sum(i);
		if(safe_isnan(mass)) {
			critical("Encountered a nan while computing the total mass (bin {})", i);
		}
	}

	return mass * V_bin / _sim_state.rho.bins();
}

template<int dims>
double CahnHilliard<dims>::average_free_energy() {
	integrator->sync();

	double fe = 0.;
	for(unsigned int i = 0; i < _sim_state.rho.bins(); i++) {
		double interfacial_contrib = 0.;
		for(int species = 0; species < model->N_species(); species++) {
			auto rho_grad = gradient(_sim_state.rho, species, i);
			for(int d = 0; d < dims; d++) {
				interfacial_contrib += k_laplacian * rho_grad[d] * rho_grad[d];
			}
		}
		fe += model->bulk_free_energy(_sim_state.rho.species_view(i)) + interfacial_contrib;
		if(safe_isnan(fe)) {
			critical("Encountered a nan while computing the total free energy (bin {})", i);
		}
	}

	return fe * V_bin / _sim_state.rho.bins();
}

template<int dims>
double CahnHilliard<dims>::average_pressure() {
	integrator->sync();

	double pressure = 0.;
	for(unsigned int i = 0; i < _sim_state.rho.bins(); i++) {
		double interfacial_contrib = 0.;
		for(int species = 0; species < model->N_species(); species++) {
			pressure += model->pressure(species, _sim_state.rho.species_view(i));
			auto rho_grad = gradient(_sim_state.rho, species, i);
			for(int d = 0; d < dims; d++) {
				pressure -= 0.5 * k_laplacian * rho_grad[d] * rho_grad[d];
			}
		}
	}

	return pressure / _sim_state.rho.bins();
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, const std::string &filename, long long int t) {
	std::ofstream output(filename);
	print_species_density(species, output, t);
	output.close();
}

template<int dims>
void CahnHilliard<dims>::print_species_density(int species, std::ofstream &output, long long int t) {
	integrator->sync();
}

template<int dims>
void CahnHilliard<dims>::print_total_density(const std::string &filename, long long int t) {
	integrator->sync();

	std::ofstream output(filename.c_str());

	output << fmt::format("# step = {}, t = {:.5}, size = {}", t, t * dt, _grid_size_str) << std::endl;
	int newline_every = (dims == 1) ? 1 : N;
	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0 && idx % newline_every == 0) {
			output << std::endl;
		}
		output << _density_to_user(_sim_state.rho.field_sum(idx)) << " ";
	}
	output << std::endl;

	output.close();
}

template<int dims>
void CahnHilliard<dims>::print_pressure(const std::string &filename, long long int t) {
	std::ofstream output(filename);
	print_pressure(output, t);
	output.close();
}

template<int dims>
void CahnHilliard<dims>::print_pressure(std::ofstream &output, long long int t) {
	integrator->sync();
	
	output << fmt::format("# pressure @ step = {} t = {:.5} size = {}", t, t * dt, _grid_size_str) << std::endl;
	int newline_every = (dims == 1) ? 1 : N;
	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0 && idx % newline_every == 0) {
			output << std::endl;
		}
		
		double pressure = 0.;
		for(int species = 0; species < model->N_species(); species++) {
			pressure += model->pressure(species, _sim_state.rho.species_view(idx));
			auto rho_grad = gradient(_sim_state.rho, species, idx);
			for(int d = 0; d < dims; d++) {
				pressure -= 0.5 * k_laplacian * rho_grad[d] * rho_grad[d];
			}
		}
		output << pressure << " ";
	}
	output << std::endl;
}

template<int dims>
double CahnHilliard<dims>::_density_to_user(double v) {
	return v / (_internal_to_user * _internal_to_user * _internal_to_user);
}

template class CahnHilliard<1>;
template class CahnHilliard<2>;

} /* namespace ch */
