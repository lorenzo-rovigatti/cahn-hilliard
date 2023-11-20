/*
 * SimpleWertheim.cuh
 *
 *  Created on: Nov 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_MODELS_SIMPLEWERTHEIM_CUH_
#define SRC_CUDA_MODELS_SIMPLEWERTHEIM_CUH_

#include "../../defs_CUDA.h"

namespace ch {
    void simple_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2, int valence, float two_valence_delta);
}

#endif /* SRC_CUDA_MODELS_SIMPLEWERTHEIM_CUH_ */