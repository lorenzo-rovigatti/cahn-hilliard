/*
 * CUDAIntegrator.h
 *
 * Created on: 3/29/2024
 *      Author: Lorenzo
*/

#ifndef CUDAINTEGRATOR_H
#define CUDAINTEGRATOR_H

#include "../../defs_CUDA.h"

#include "../../integrators/Integrator.h"

namespace ch {

template<int dims>
class CUDAIntegrator : public Integrator<dims> {
public:
    CUDAIntegrator(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config);

    ~CUDAIntegrator();

    void sync() override;

    GET_NAME(CUDAIntegrator)

protected:
    void _CPU_GPU();

    int _d_vec_size;
	MultiField<field_type> _h_rho;
	field_type *_d_rho = nullptr; // points to sim_state.CUDA_rho
	float *_d_rho_der = nullptr;

    MultiField<field_type> _h_mobility;
	field_type *_d_mobility = nullptr; // points to sim_state.CUDA_mobility

    int _grid_size;
    bool _output_ready = false;
};

} /* namespace ch */

#endif /* CUDA_INTEGRATOR_H */
