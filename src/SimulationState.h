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

template<int dims> class Integrator;

template<int dims>
struct SimulationState {
    MultiField<double> rho;
    MultiField<double> mobility;

#ifndef NOCUDA
    field_type *CUDA_mobility = nullptr;
    field_type *CUDA_rho = nullptr;
#endif

    std::unique_ptr<ch::FreeEnergyModel> model;
    std::unique_ptr<ch::Integrator<dims>> *integrator;

    double user_to_internal;
    double internal_to_user;

    bool use_CUDA;

    double density_to_user(double v) const {
        return v / (internal_to_user * internal_to_user * internal_to_user);
    }

    double user_to_density(double v) const {
        return v * (user_to_internal * user_to_internal * user_to_internal);
    }

    double length_to_user(double v) const {
        return v * internal_to_user;
    }

    double user_to_length(double v) const {
        return v * user_to_internal;
    }
};

} // namespace ch

#endif /* SIMULATIONSTATE_H */
