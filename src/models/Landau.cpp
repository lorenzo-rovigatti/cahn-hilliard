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

	
	if(config["landau"]["epsilon"]) {
		_epsilon = _config_value<double>(config, "landau.epsilon");
	}
	else if(config["landau"]["T"]) {
		double T = _config_value<double>(config, "landau.T");
		_epsilon = 1 - T;
	}
	else {
		critical("The Landau free energy requires either 'epsilon' or 'T' to be set");
	}

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

void Landau::der_bulk_free_energy(const MultiField<double> &rho, MultiField<double> &rho_der) {
	for(unsigned int idx = 0; idx < rho.bins(); idx++) {
        for(int species = 0; species < N_species(); species++) {
			double op = rho.species_view(idx)[0];
            rho_der(idx, species) = -_epsilon * op + op * op * op;
        }
    }
}

double Landau::bulk_free_energy(const SpeciesView<double> &rhos) {
	double op = rhos[0];
	return -0.5 * _epsilon * SQR(op) + 0.25 * SQR(SQR(op));
}

void Landau::der_bulk_free_energy(field_type *psi, float *psi_der, int grid_size) {
#ifndef NOCUDA
	landau_der_bulk_free_energy(psi, psi_der, grid_size, _epsilon);
#endif
}

double Landau::pressure(int species, const SpeciesView<double> &rhos) {
	double op = rhos[species];
	double mu = -_epsilon * op + op * op * op;
	double f = bulk_free_energy(rhos);
	return op * mu - f;
}

} /* namespace ch */
