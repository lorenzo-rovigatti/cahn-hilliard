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
    CUDAIntegrator(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);

    ~CUDAIntegrator();

    MultiField<double> &rho();

    GET_NAME(CUDAIntegrator)

protected:
    void _CPU_GPU();

    int _d_vec_size;
	MultiField<field_type> _h_rho;
	field_type *_d_rho = nullptr;
	float *_d_rho_der = nullptr;
    int _grid_size;
    bool _output_ready = false;
};

} /* namespace ch */

#endif /* CUDA_INTEGRATOR_H */
