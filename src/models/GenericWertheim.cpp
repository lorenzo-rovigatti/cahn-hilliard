/*
 * GenericWertheim.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: lorenzo
 */

#include "GenericWertheim.h"
#include "../utils/Delta.h"

#ifndef NOCUDA
// #include "../CUDA/models/GenericWertheim.cuh"
#endif

#include <numeric>

namespace ch {

GenericWertheim::GenericWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_B2 = _config_value<double>(config, "generic_wertheim.B2");
	_B3 = _config_optional_value<double>(config, "generic_wertheim.B3", 0.);
	_delta_AA = Delta(config, "generic_wertheim.delta_AA");
	_delta_BB = Delta(config, "generic_wertheim.delta_BB");
	_delta_CC = Delta(config, "generic_wertheim.delta_CC");
	// _valence = _config_array_values<int>(config, "generic_wertheim.valence", 3);
	_valence = _config_value<int>(config, "generic_wertheim.valence");
	_linker_partial_valence = 2;

	info("valence = {}, linker_partial_valence = {}, B2 = {}, delta_AA = {}, delta_BB = {}", _valence, _linker_partial_valence, _B2, _delta_AA, _delta_BB);

	_B2 *= CUB(_user_to_internal);
	_B3 *= SQR(CUB(_user_to_internal));
	_delta_AA *= CUB(_user_to_internal);
	_delta_BB *= CUB(_user_to_internal);

#ifndef NOCUDA
	// if(_config_optional_value<bool>(config, "use_CUDA", false)) {
	// 	init_Generic_symbols(_valence, _linker_half_valence, _delta_AA, _delta_BB);
	// }
#endif
}

GenericWertheim::~GenericWertheim() {
}

double GenericWertheim::bonding_free_energy(const std::vector<double> &rhos) {
	double rho_factor =  _delta_AA * (_valence * rhos[0] + _linker_partial_valence * rhos[3]);
	double X_1A = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_1 = std::log(X_1A) - X_1A / 2.0 + 0.5;

	rho_factor =  _delta_BB * (_valence * rhos[1] + _linker_partial_valence * rhos[3]);
	double X_2B = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_2 = std::log(X_2B) - X_2B / 2.0 + 0.5;

	rho_factor =  _delta_CC * (_valence * rhos[2] + _linker_partial_valence * rhos[3]);
	double X_3C = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_3 = std::log(X_3C) - X_3C / 2.0 + 0.5;

	double bonding_fe =
			rhos[0] * _valence * fe_part_1 +
			rhos[1] * _valence * fe_part_2 +
			rhos[2] * _valence * fe_part_3 +
			rhos[3] * _linker_partial_valence * (fe_part_1 + fe_part_2 + fe_part_3);

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

double GenericWertheim::_der_contribution(const std::vector<double> &rhos, int species) {
	double delta = (species == 0) ? _delta_AA : _delta_BB;
	double rho_factor =  delta * (_valence * rhos[species] + _linker_partial_valence * rhos[3]);
	double X = (-1.0 + std::sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	return (rho_factor >= 0) ? std::log(X) : 0.0;
}

double GenericWertheim::der_bulk_free_energy(int species, const std::vector<double> &rhos) {
	// early return if the species has zero density
	if(rhos[species] == 0.) {
		return 0.;
	}
	
	double rho = std::accumulate(rhos.begin(), rhos.end(), 0.);
	double der_f_ref = std::log(rhos[species]) + 2.0 * _B2 * rho + 3.0 * _B3 * SQR(rho);

	double der_f_bond;
	if(species == 3) {
		der_f_bond = _linker_partial_valence * (_der_contribution(rhos, 0) + _der_contribution(rhos, 1) + _der_contribution(rhos, 2));
	}
	else {
		der_f_bond = _valence * _der_contribution(rhos, species);
	}

	return der_f_ref + der_f_bond;
}

void GenericWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
#ifndef NOCUDA
	// Generic_wertheim_der_bulk_free_energy(rho, rho_der, vec_size, _B2, _B3);
#endif
}

} /* namespace ch */
