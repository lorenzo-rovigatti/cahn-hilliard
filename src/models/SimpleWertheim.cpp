/*
 * SimpleWertheim.cpp
 *
 *  Created on: Jul 22, 2023
 *      Author: lorenzo
 */

#include "SimpleWertheim.h"

#include "../utils/utility_functions.h"
#include "../utils/Delta.h"

#ifndef NOCUDA
#include "../CUDA/models/SimpleWertheim.cuh"
#endif

#include <iostream>

namespace ch {

SimpleWertheim::SimpleWertheim(toml::table &config) : FreeEnergyModel(config) {
    _valence = _config_value<int>(config, "wertheim.valence");
    _B2 = _config_value<double>(config, "wertheim.B2");
    _regularisation_delta = _config_optional_value<double>(config, "wertheim.regularisation_delta", 0.0);
    _delta = Delta(config, "wertheim.delta");

    info("valence = {}, delta = {}", _valence, _delta);

    _B2 *= CUB(_user_to_internal);
    _delta *= CUB(_user_to_internal);
    _regularisation_delta /= CUB(_user_to_internal);

    _log_delta = std::log(_regularisation_delta);
    _two_valence_delta = 2.0 * _valence * _delta;
}

SimpleWertheim::~SimpleWertheim() {

}

double SimpleWertheim::der_bulk_free_energy_expansive(int species, const std::vector<double> &rhos) {
    double rho = rhos[species];
    double X = _X(rho);
    double der_f_bond = (rho > 0.) ? _valence * std::log(X) : 0.0;

    return der_f_bond;
}

double SimpleWertheim::der_bulk_free_energy_contractive(int species, const std::vector<double> &rhos) {
    double rho = rhos[species];
    double der_f_ref = (rho < _regularisation_delta) ? rho / _regularisation_delta + _log_delta - 1.0 : std::log(rho);
	der_f_ref += 2 * _B2 * rho;

    return der_f_ref;
}

double SimpleWertheim::der_bulk_free_energy(int species, const std::vector<double> &rhos) {
    double rho = rhos[species];
    double der_f_ref = (rho < _regularisation_delta) ? rho / _regularisation_delta + _log_delta - 1.0 : std::log(rho);
	der_f_ref += 2 * _B2 * rho;
    double X = _X(rho);
    double der_f_bond = (rho > 0.) ? _valence * std::log(X) : 0.0;

    return (der_f_ref + der_f_bond);
}

double SimpleWertheim::bulk_free_energy(const std::vector<double> &rhos) {
    double rho = rhos[0];
    double f_ref = (rho < _regularisation_delta) ? SQR(rho) / (2.0 * _regularisation_delta) + rho * _log_delta - _regularisation_delta / 2.0 : rho * std::log(rho * _density_conversion_factor);
    f_ref += -rho + _B2 * SQR(rho);
    double f_bond = (rho > 0.) ? _valence * rho * (std::log(_X(rho)) + 0.5 * (1. - _X(rho))) : 0.0;

    return (f_ref + f_bond);
}

double SimpleWertheim::_X(double rho) {
    double sqrt_argument = 2.0 * _two_valence_delta * rho;
    return (-1.0 + std::sqrt(1.0 + 2.0 * _two_valence_delta * rho)) / (_two_valence_delta * rho);
}

void SimpleWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size) {
#ifndef NOCUDA
    simple_wertheim_der_bulk_free_energy(rho, rho_der, grid_size, _B2, _valence, _two_valence_delta);
#endif
}

} /* namespace ch */
