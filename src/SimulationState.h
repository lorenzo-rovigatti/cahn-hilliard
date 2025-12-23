/*
 * SimulationState.h
 *
 * Created on: 12/19/2025
 *      Author: Lorenzo
*/

#ifndef SIMULATIONSTATE_H
#define SIMULATIONSTATE_H

#include "utils/MultiField.h"
#include "models/FreeEnergyModel.h"

namespace ch {

struct SimulationState {
    MultiField<double> rho;
    MultiField<double> mobility;

    std::unique_ptr<ch::FreeEnergyModel> model;

    double user_to_internal;
    double internal_to_user;

    long long int time_step;
};

} // namespace ch

#endif /* SIMULATIONSTATE_H */
