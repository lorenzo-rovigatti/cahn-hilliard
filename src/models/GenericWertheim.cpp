/*
 * GenericWertheim.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: lorenzo
 */

#include "GenericWertheim.h"
#include "../utils/Delta.h"
#include "../utils/strings.h"

#ifndef NOCUDA
// #include "../CUDA/models/GenericWertheim.cuh"
#endif

#include <numeric>

namespace ch {

GenericWertheim::GenericWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_B2 = _config_value<double>(config, "generic_wertheim.B2");
	_B3 = _config_optional_value<double>(config, "generic_wertheim.B3", 0.);

	_B2 *= CUB(_user_to_internal);
	_B3 *= SQR(CUB(_user_to_internal));

	if(auto arr = config["generic_wertheim"]["species"].as_array()) {
		std::set<int> unique_patches_set;
        for(const auto& elem : *arr) {
			Species new_species;
            new_species.patches = _config_array_values<int>(*elem.as_table(), "patches");
			new_species.idx = _species.size();

			// here we fill the two unique_* vectors by using an auxiliary map
			std::map<int, int> counter;
			for(int val : new_species.patches) {
				counter[val]++;
			}

			for (const auto& [key, count] : counter) {
				UniquePatch up {
					.idx = key, 
					.multiplicity = count};
				new_species.unique_patches.push_back(up);
			}
			new_species.N_unique_patches = new_species.unique_patches.size();

			_species.push_back(new_species);
			unique_patches_set.insert(new_species.patches.begin(), new_species.patches.end());
        }

		_unique_patch_ids.assign(unique_patches_set.begin(), unique_patches_set.end());
		_N_patches = *std::max_element(_unique_patch_ids.begin(), _unique_patch_ids.end()) + 1;
		_delta.resize(_N_patches * _N_patches, 0.0);
    }
	else {
        critical("Missing 'generic_wertheim.species' array");
    }

	if(auto arr = config["generic_wertheim"]["deltas"].as_array()) {
        for(const auto& elem : *arr) {
			auto interaction = _config_value<std::string>(*elem.as_table(), "interaction");
			auto split = utils::split(interaction, "-");
			if(split.size() != 2) {
				critical("The following delta interaction specifier is malformed: {}", interaction);
			}
			// TODO: check that these are effectively numbers
			int patch_A = std::atoi(split[0].c_str());
			int patch_B = std::atoi(split[1].c_str());
			double my_delta = Delta(toml::node_view<const toml::node>{ elem }); // the horror
			my_delta *= CUB(_user_to_internal);
			_delta[patch_A * _N_patches + patch_B] = _delta[patch_B * _N_patches + patch_A] = my_delta;
			info("Adding {}-{} interaction with delta = {:.2f}", patch_A, patch_B, my_delta);
		}
    }
	else {
        critical("Missing 'generic_wertheim.deltas' array");
    }

	// build the list of patch interactions for each unique patch
	_unique_patch_interactions.resize(_N_patches);
	for(size_t i = 0; i < _unique_patch_ids.size(); i++) {
		int patch = _unique_patch_ids[i];
		std::vector<PatchInteraction> patch_interacting_species;
		for(const auto& other : _species) {
			PatchInteraction interaction;
			for(auto &other_patch : other.unique_patches) {
				if(_delta[patch * _N_patches +  other_patch.idx] > 0.0) {
					interaction.species = other.idx;
					interaction.patches.push_back(other_patch);
				}
			}
			if(interaction.species != -1) {
				patch_interacting_species.push_back(interaction);
			}
		}
		_unique_patch_interactions[patch] = patch_interacting_species;
	}

#ifndef NOCUDA
	// if(_config_optional_value<bool>(config, "use_CUDA", false)) {
	// 	init_Generic_symbols(_valence, _linker_half_valence, _delta_AA, _delta_BB);
	// }
#endif
}

GenericWertheim::~GenericWertheim() {
}

void GenericWertheim::_update_X(const std::vector<double> &rhos, std::vector<double> &Xs) {
	double tolerance = 1e-8;
	int max_iter = 10000;

	for(int iter = 0; iter < max_iter; iter++) {
		double max_delta = 0.0;

		for(auto &patch : _unique_patch_ids) {
			double sum = 0.0;
			for(auto &interaction : _unique_patch_interactions[patch]) {
				double rho = rhos[interaction.species];
				for(auto &other_patch : interaction.patches) {
					double delta = _delta[patch * _N_patches +  other_patch.idx];
					sum += other_patch.multiplicity * rho * Xs[other_patch.idx] * delta;
				}
				
			}

			double new_X = 1.0 / (1.0 + sum);
			max_delta = std::max(max_delta, std::fabs(new_X - Xs[patch]));
			Xs[patch] = new_X;
		}

		if(max_delta < tolerance) {
			return;
		};
	}
}

double GenericWertheim::bonding_free_energy(const std::vector<double> &rhos) {
	double bonding_fe = 0;
	std::vector<double> Xs(_N_patches, 0.0);
	_update_X(rhos, Xs);

	for(auto &species : _species) {
		double species_fe = 0;
		for(auto &patch : species.unique_patches) {
			species_fe += patch.multiplicity * (std::log(Xs[patch.idx]) - Xs[patch.idx] / 2.0 + 0.5);
		}
		bonding_fe += rhos[species.idx] * species_fe;
	}
	
	return bonding_fe;
}

double GenericWertheim::bulk_free_energy(const std::vector<double> &rhos) {
	double rho = std::accumulate(rhos.begin(), rhos.end(), 0.);

	double mixing_S = 0.;
	for(int i = 0; i < N_species(); i++) {
		double x_i = rhos[i] / rho;
		mixing_S += (x_i > 0) ? x_i * std::log(x_i) : 0;
	}
	double B2_contrib = _B2 * rho;
	double B3_contrib = 0.5 * _B3 * SQR(rho);

	double f_ref = rho * (std::log(rho * _density_conversion_factor) - 1.0 + mixing_S + B2_contrib + B3_contrib);

	return f_ref + bonding_free_energy(rhos);
}

void GenericWertheim::der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) {
	static std::vector<std::vector<double>> all_Xs(rho.bins(), std::vector<double>(_N_patches, 0.0));

	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		std::vector<double> rhos = rho.rho_species(idx);
		std::vector<double> &Xs = all_Xs[idx];
		_update_X(rhos, Xs);
        for(int species = 0; species < N_species(); species++) {
			// early return if the species has zero density
			if(rhos[species] != 0.) {
				double rho_tot = std::accumulate(rhos.begin(), rhos.end(), 0.);
				double der_f_ref = std::log(rhos[species]) + 2.0 * _B2 * rho_tot + 3.0 * _B3 * SQR(rho_tot);
				
				double der_f_bond = 0;
				for(auto &patch : _species[species].unique_patches) {
					der_f_bond += patch.multiplicity * std::log(Xs[patch.idx]);
				}
				rho_der(idx, species) = der_f_ref + der_f_bond;
			}
        }
    }
}

void GenericWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
#ifndef NOCUDA
	// Generic_wertheim_der_bulk_free_energy(rho, rho_der, vec_size, _B2, _B3);
#endif
}

} /* namespace ch */
