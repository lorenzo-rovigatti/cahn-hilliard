/*
 * PseudospectralCUDA.h
 *
 * Created on: 3/29/2024
 *      Author: Lorenzo
*/

#ifndef PSEUDOSPECTRALCUDA_H
#define PSEUDOSPECTRALCUDA_H

#include "../../defs_CUDA.h"

#include "CUDAIntegrator.h"

namespace ch {

template<int dims>
class PseudospectralCUDA : public CUDAIntegrator<dims> {
public:
    PseudospectralCUDA(FreeEnergyModel *model, toml::table &config);

    ~PseudospectralCUDA();

    void set_initial_rho(MultiField<double> &r) override;

    void evolve() override;

    GET_NAME(PseudospectralCUDA)

private:
    int _hat_vector_size;
    cufftFieldComplex *_d_rho_hat = nullptr, *_d_rho_hat_copy;
	cufftComplex *_d_f_der_hat = nullptr;
	float *_d_sqr_wave_vectors = nullptr; 
	float *_d_dealiaser = nullptr;

	cufftHandle _d_rho_inverse_plan, _d_f_der_plan;
};

} /* namespace ch */

#endif /* PSEUDOSPECTRALCUDA_H */
