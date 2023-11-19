/*
 * Landau.cpp
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#include "Landau.h"

#include "../CUDA/models/Landau.cuh"

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

double Landau::bulk_free_energy(std::vector<double> &rhos) {
	double op = rhos[0];
	return -0.5 * _epsilon * SQR(op) + 0.25 * SQR(SQR(op));
}

void Landau::der_bulk_free_energy(double *psi, float *psi_der, int grid_size) {
	landau_der_bulk_free_energy(psi, psi_der, grid_size, _epsilon);
}

} /* namespace ch */
