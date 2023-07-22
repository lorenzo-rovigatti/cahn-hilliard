/*
 * SimpleWertheim.cpp
 *
 *  Created on: Jul 22, 2023
 *      Author: lorenzo
 */

#include "SimpleWertheim.h"

namespace ch {

SimpleWertheim::SimpleWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_valence = _config_value<int>(config, "wertheim.valence");
	_B2 = _config_value<double>(config, "wertheim.B2");
	_regularisation_delta = _config_optional_value<double>(config, "wertheim.regularisation_delta", 0.0);
	_log_delta = std::log(_regularisation_delta);
	_T = _config_value<double>(config, "wertheim.T");

	double salt = _config_optional_value<double>(config, "wertheim.salt", 1.0);
	int L_DNA = _config_optional_value<int>(config, "wertheim.sticky_size", 6);
	double delta_H = _config_value<double>(config, "wertheim.deltaH");
	double delta_S = _config_value<double>(config, "wertheim.deltaS");

	double delta_S_salt = 0.368 * (L_DNA - 1.0) * std::log(salt);
	double delta_G = delta_H - _T * (delta_S + delta_S_salt);

	const double k_B = 1.9872036;
	_delta = 1.6606 * std::exp(-delta_G / (k_B * _T));
	_two_valence_delta = 2 * _valence * _delta;

	info("valence = {}, delta = {}", _valence, _delta);
}

SimpleWertheim::~SimpleWertheim() {

}

double SimpleWertheim::bulk_free_energy(int species, std::vector<double> &rhos) {
	double rho = rhos[species];
	double f_ref = rho * std::log(rho) - rho + _B2 * SQR(rho);
	double f_bond = _valence * rho * (std::log(_X(rho)) + 0.5 * (1. - _X(rho)));

	return (f_ref + f_bond);
}

double SimpleWertheim::der_bulk_free_energy(int species, std::vector<double> &rhos) {
	double rho = rhos[species];
	double der_f_ref = _regularised_log(rho) + 2 * _B2 * rho;
	double X = _X(rho);
	double der_f_bond = _valence * (std::log(X) - 0.5 + 1.0 / (2.0 - X) - 0.5 * X / (2.0 - X));

	return (der_f_ref + der_f_bond);
}

double SimpleWertheim::_X(double rho) {
	double sqrt_argument = 2.0 * _two_valence_delta * rho;
	return (sqrt_argument < 1e-3) ? 1.0 - sqrt_argument / 2.0 : (-1 + std::sqrt(1 + 2 * _two_valence_delta * rho)) / (_two_valence_delta * rho);
}

} /* namespace ch */
