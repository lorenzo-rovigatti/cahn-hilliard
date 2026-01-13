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

#ifndef NOCUDA
    field_type *CUDA_mobility = nullptr;
    field_type *CUDA_rho = nullptr;
#endif

    std::unique_ptr<ch::FreeEnergyModel> model;

    double user_to_internal;
    double internal_to_user;

    long long int time_step;

    bool use_CUDA;
};

} // namespace ch

#endif /* SIMULATIONSTATE_H */
