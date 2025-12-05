/*
 * SalehWertheim.cpp
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#include "SalehWertheim.h"
#include "../utils/Delta.h"

#ifndef NOCUDA
#include "../CUDA/models/SalehWertheim.cuh"
#endif

#include <numeric>

namespace ch {

SalehWertheim::SalehWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_B2 = _config_value<double>(config, "saleh.B2");
	_B3 = _config_optional_value<double>(config, "saleh.B3", 0.);
	_delta_AA = Delta(config, "saleh.delta_AA");
	_delta_BB = Delta(config, "saleh.delta_BB");
	_valence = _config_array_values<int>(config, "saleh.valence", 3);
	_linker_half_valence = _valence[2] / 2;

	info("valences = ({}), B2 = {}, delta_AA = {}, delta_BB = {}", fmt::join(_valence, ", "), _B2, _delta_AA, _delta_BB);

	_B2 *= CUB(_user_to_internal);
	_B3 *= SQR(CUB(_user_to_internal));
	_delta_AA *= CUB(_user_to_internal);
	_delta_BB *= CUB(_user_to_internal);

#ifndef NOCUDA
	if(_config_optional_value<bool>(config, "use_CUDA", false)) {
		init_saleh_symbols(_valence, _linker_half_valence, _delta_AA, _delta_BB);
	}
#endif
}

SalehWertheim::~SalehWertheim() {
}

double SalehWertheim::bonding_free_energy(const SpeciesView<double> &rhos) {
	double rho_factor =  _delta_AA * (_valence[0] * rhos[0] + _linker_half_valence * rhos[2]);
	double X_1A = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_1 = std::log(X_1A) - X_1A / 2.0 + 0.5;

	rho_factor =  _delta_BB * (_valence[1] * rhos[1] + _linker_half_valence * rhos[2]);
	double X_2B = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_2 = std::log(X_2B) - X_2B / 2.0 + 0.5;

	double bonding_fe =
			rhos[0] * _valence[0] * fe_part_1 +
			rhos[1] * _valence[1] * fe_part_2 +
			rhos[2] * _linker_half_valence * (fe_part_1 + fe_part_2);

	return bonding_fe;
}

double SalehWertheim::bulk_free_energy(const SpeciesView<double> &rhos) {
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

double SalehWertheim::_der_contribution(const SpeciesView<double> &rhos, int species) {
	double delta = (species == 0) ? _delta_AA : _delta_BB;
	double rho_factor =  delta * (_valence[species] * rhos[species] + _linker_half_valence * rhos[2]);
	double X = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	return (rho_factor >= 0) ? std::log(X) : 0.0;
}

void SalehWertheim::der_bulk_free_energy(const MultiField<double> &rho, MultiField<double> &rho_der) {
	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
        for(int species = 0; species < N_species(); species++) {
			auto rhos = rho.species_view(idx);
			if(rhos[species] != 0.) {
				double rho_tot = std::accumulate(rhos.begin(), rhos.end(), 0.);
				double der_f_ref = std::log(rhos[species]) + 2.0 * _B2 * rho_tot + 3.0 * _B3 * SQR(rho_tot);

				double der_f_bond;
				if(species == 0) {
					der_f_bond = _valence[0] * _der_contribution(rhos, species);
				}
				else if(species == 1) {
					der_f_bond = _valence[1] * _der_contribution(rhos, species);
				}
				else {
					der_f_bond = _linker_half_valence * (_der_contribution(rhos, 0) + _der_contribution(rhos, 1));
				}
				rho_der(idx, species) = der_f_ref + der_f_bond;
			}
        }
    }
}

void SalehWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
#ifndef NOCUDA
	saleh_wertheim_der_bulk_free_energy(rho, rho_der, vec_size, _B2, _B3);
#endif
}

} /* namespace ch */
