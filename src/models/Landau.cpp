/*
 * Landau.cpp
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#include "../models/Landau.h"

namespace ch {

Landau::Landau(toml::table &config) :
				FreeEnergyModel(config) {

	_epsilon = _config_value<double>(config, "landau.epsilon");

	info("epsilon = {}", _epsilon);
}

Landau::~Landau() {

}

int Landau::N_species() {
	return 1;
}

double Landau::der_bulk_free_energy(int species, std::vector<double> &rhos) {
	double op = rhos[species];
	return -_epsilon * op + op * op * op;
}

double Landau::bulk_free_energy(int species, std::vector<double> &rhos) {
	double op = rhos[species];
	return -0.5 * _epsilon * SQR(op) + 0.25 * SQR(SQR(op));
}

} /* namespace ch */
