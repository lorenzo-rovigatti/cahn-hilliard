/*
 * SalehWertheim.cuh
 *
 *  Created on: Nov 21, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_MODELS_SALEHWERTHEIM_CUH_
#define SRC_CUDA_MODELS_SALEHWERTHEIM_CUH_

#include <vector>

#include "../../defs_CUDA.h"

namespace ch {
    void init_saleh_symbols(std::vector<int> &valence, int linker_half_valence, float delta_AA, float delta_BB);
    void saleh_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2);
}

#endif /* SRC_CUDA_MODELS_SALEHWERTHEIM_CUH_ */