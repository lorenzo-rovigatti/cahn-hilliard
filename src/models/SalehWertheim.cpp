/*
 * SalehWertheim.cpp
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#include "SalehWertheim.h"
#include "../utils/Delta.h"

#include <numeric>

namespace ch {

SalehWertheim::SalehWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_B2 = _config_value<double>(config, "saleh.B2");
	_delta_AA = Delta(config, "saleh.delta_AA");
	_delta_BB = Delta(config, "saleh.delta_BB");
	_valence = _config_array_values<int>(config, "saleh.valence", 3);
	_linker_half_valence = _valence[2] / 2;

	info("valences = ({}), B2 = {}, delta_AA = {}, delta_BB = {}", fmt::join(_valence, ", "), _B2, _delta_AA, _delta_BB);
}

SalehWertheim::~SalehWertheim() {
}

double SalehWertheim::bonding_free_energy(std::vector<double> &rhos) {
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

double SalehWertheim::bulk_free_energy(std::vector<double> &rhos) {
	double rho = std::accumulate(rhos.begin(), rhos.end(), 0.);

	double mixing_S = 0., B2_contrib = 0.;
	for(int i = 0; i < N_species(); i++) {
		double x_i = rhos[i] / rho;
		mixing_S += x_i * std::log(x_i);

		for(int j = 0; j < N_species(); j++) {
			B2_contrib += _B2 * x_i * rhos[j];
		}
	}

	double f_ref = rho * (std::log(rho) - 1.0 + mixing_S + B2_contrib);

	return f_ref + bonding_free_energy(rhos);
}

double SalehWertheim::der_bulk_free_energy(int species, std::vector<double> &rhos) {
	// the ideal + B2 part is computed analytically
	double der_f_ref = std::log(rhos[species]);
	for(int i = 0; i < N_species(); i++) {
		der_f_ref += 2.0 * _B2 * rhos[i];
	}

	// the bonding part is computed numerically
	double rho_i = rhos[species];
	double delta_rho_i = rho_i * 1e-5;

	double fe_r = bonding_free_energy(rhos);
	rhos[species] += delta_rho_i;
	double fe_rdr = bonding_free_energy(rhos);
	rhos[species] = rho_i;

	return der_f_ref + (fe_rdr - fe_r) / delta_rho_i;
}

} /* namespace ch */