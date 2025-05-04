/*
 * RicciWertheim.cuh
 *
 *  Created on: May 4, 2025
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_MODELS_RICCIWERTHEIM_CUH_
#define SRC_CUDA_MODELS_RICCIWERTHEIM_CUH_

#include <vector>

#include "../../defs_CUDA.h"

namespace ch {
    void init_ricci_symbols(float delta_00, float delta_12);
    void ricci_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2);
}

#endif /* SRC_CUDA_MODELS_RICCIWERTHEIM_CUH_ */