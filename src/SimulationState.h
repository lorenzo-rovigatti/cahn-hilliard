/*
 * SimulationState.h
 *
 * Created on: 12/19/2025
 *      Author: Lorenzo
*/

#ifndef SIMULATIONSTATE_H
#define SIMULATIONSTATE_H

#include "utils/MultiField.h"

namespace ch {

struct SimulationState {
    MultiField<double> rho;
    MultiField<double> mobility;

    long long int time_step;
};

} // namespace ch

#endif /* SIMULATIONSTATE_H */
