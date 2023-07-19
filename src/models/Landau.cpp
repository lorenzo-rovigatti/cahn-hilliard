/*
 * Landau.cpp
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#include "../models/Landau.h"

namespace ch {

Landau::Landau(cxxopts::Options &options) :
				FreeEnergyModel(options) {
	options.add_options()
	("e,epsilon", "The distance from the critical point in the 'landau' free energy", cxxopts::value<double>()->default_value("0.9"));
}

Landau::~Landau() {

}

void Landau::init(cxxopts::ParseResult &result) {
	_epsilon = result["epsilon"].as<double>();
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
