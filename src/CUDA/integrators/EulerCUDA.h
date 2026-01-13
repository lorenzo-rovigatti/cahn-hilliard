/*
 * EulerCUDA.h
 *
 * Created on: 3/29/2024
 *      Author: Lorenzo
*/

#ifndef EULERCUDA_H
#define EULERCUDA_H

#include "../../defs_CUDA.h"

#include "CUDAIntegrator.h"

namespace ch {

template<int dims>
class EulerCUDA : public CUDAIntegrator<dims> {
public:
    EulerCUDA(SimulationState &sim_state,FreeEnergyModel *model, toml::table &config);

    ~EulerCUDA();

    void evolve() override;

    GET_NAME(EulerCUDA)
};

} /* namespace ch */

#endif /* EULERCUDA_H */
