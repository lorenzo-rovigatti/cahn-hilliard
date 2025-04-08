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
        for(const auto& elem : *arr) {
			Species new_species;
            new_species.patches = _config_array_values<int>(*elem.as_table(), "patches");
			new_species.N_patches = new_species.patches.size();
			new_species.idx = _species.size();
			_species.push_back(new_species);
        }

    } else {
        critical("Key 'generic_wertheim.species' should be an array\n");
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
			// TODO: add log printing of the parsed values
			// _delta[{patch_A, patch_B}] = Delta(*elem.as_table(), );
		}
    } else {
        critical("Key 'generic_wertheim.deltas' should be an array\n");
    }

	_delta[{0, 0}] = 10000;

	for(auto &species : _species) {
		_N_patches += species.N_patches;
	}

	for(auto &delta : _delta) {
		delta.second *= CUB(_user_to_internal);
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
	int max_iter = 1000;

	for(int iter = 0; iter < max_iter; ++iter) {
		double max_delta = 0.0;

		for(auto& species : _species) {
			for(auto &patch : species.patches) {
				double sum = 0.0;
				for(const auto& other : _species) {
					for(auto &other_patch : other.patches) {
						double delta = _delta[{patch, other_patch}];
						sum += other.idx * Xs[other_patch] * delta;
					}
				}

				double new_X = 1.0 / (1.0 + sum);
				max_delta = std::max(max_delta, std::abs(new_X - Xs[patch]));
				Xs[patch] = new_X;
			}
		}

		if(max_delta < tolerance) break;
	}
}

double GenericWertheim::bonding_free_energy(const std::vector<double> &rhos) {
	double bonding_fe = 0;
	std::vector<double> Xs(_N_patches, 0.0);
	_update_X(rhos, Xs);

	for(auto &species : _species) {
		double species_fe = 0;
		for(auto &patch : species.patches) {
			species_fe += std::log(Xs[patch]) - Xs[patch] / 2.0;
		}
		species_fe += species.N_patches / 2.0;
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

double GenericWertheim::der_bulk_free_energy(int species, const std::vector<double> &rhos) {
	// early return if the species has zero density
	if(rhos[species] == 0.) {
		return 0.;
	}
	
	double rho = std::accumulate(rhos.begin(), rhos.end(), 0.);
	double der_f_ref = std::log(rhos[species]) + 2.0 * _B2 * rho + 3.0 * _B3 * SQR(rho);

	std::vector<double> Xs(_N_patches, 0.0);
	_update_X(rhos, Xs);
	double der_f_bond = 0;
	for(auto &patch : _species[species].patches) {
		der_f_bond += std::log(Xs[patch]);
	}

	return der_f_ref + der_f_bond;
}

void GenericWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
#ifndef NOCUDA
	// Generic_wertheim_der_bulk_free_energy(rho, rho_der, vec_size, _B2, _B3);
#endif
}

} /* namespace ch */
