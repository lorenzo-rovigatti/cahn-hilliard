/*
 * Landau.cpp
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#include "Landau.h"

#ifndef NOCUDA
#include "../CUDA/models/Landau.cuh"
#endif

namespace ch {

Landau::Landau(toml::table &config) :
				FreeEnergyModel(config) {

	_epsilon = _config_value<double>(config, "landau.epsilon");

	if(_user_to_internal != 1.0) {
		critical("The Landau free energy does not support distance_scaling_factor values different from 1.0");
	}

	info("epsilon = {}", _epsilon);
}

Landau::~Landau() {

}

int Landau::N_species() {
	return 1;
}

double Landau::der_bulk_free_energy(int species, const std::vector<double> &rhos) {
	double op = rhos[species];
	return -_epsilon * op + op * op * op;
}

double Landau::bulk_free_energy(const std::vector<double> &rhos) {
	double op = rhos[0];
	return -0.5 * _epsilon * SQR(op) + 0.25 * SQR(SQR(op));
}

void Landau::der_bulk_free_energy(field_type *psi, float *psi_der, int grid_size) {
#ifndef NOCUDA
	landau_der_bulk_free_energy(psi, psi_der, grid_size, _epsilon);
#endif
}

} /* namespace ch */
