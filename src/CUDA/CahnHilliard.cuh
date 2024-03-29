/*
 * CahnHilliard.cuh
 *
 *  Created on: Nov 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_CAHNHILLIARD_CUH_
#define SRC_CUDA_CAHNHILLIARD_CUH_

#include "../defs_CUDA.h"

namespace ch {
    void init_symbols(int N, int size, int N_species);
    template<int dims> void add_surface_term(field_type *rho, float *rho_der, float dx, float k_laplacian);
    template<int dims> void integrate(field_type *rho, float *rho_der, float dx, float dt, float M);
    void integrate_fft(cufftFieldComplex *rho_hat, cufftFieldComplex *rho_hat_for_inverse_transform, cufftComplex *f_der_hat, float *sqr_wave_vectors, float *dealiaser, float dt, float M, float k_laplacian, int hat_size);
}

#endif /* SRC_CUDA_CAHNHILLIARD_CUH_ */