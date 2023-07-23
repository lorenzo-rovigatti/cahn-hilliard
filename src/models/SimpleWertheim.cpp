/*
 * SimpleWertheim.cpp
 *
 *  Created on: Jul 22, 2023
 *      Author: lorenzo
 */

#include "SimpleWertheim.h"

#include "../utils/Delta.h"

#include <iostream>

namespace ch {

SimpleWertheim::SimpleWertheim(toml::table &config) :
				FreeEnergyModel(config) {

	_valence = _config_value<int>(config, "wertheim.valence");
	_B2 = _config_value<double>(config, "wertheim.B2");
	_regularisation_delta = _config_optional_value<double>(config, "wertheim.regularisation_delta", 0.0);
	_delta = Delta(config, "wertheim.delta");

	_log_delta = std::log(_regularisation_delta);
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
