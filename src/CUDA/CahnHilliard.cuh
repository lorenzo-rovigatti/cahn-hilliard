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
    template<int dims> void add_surface_term(number *rho, number *rho_der, number dx, number k_laplacian);
    template<int dims> void integrate(number *rho, number *rho_der, number dx, number dt, number M);
}

#endif /* SRC_CUDA_CAHNHILLIARD_CUH_ */