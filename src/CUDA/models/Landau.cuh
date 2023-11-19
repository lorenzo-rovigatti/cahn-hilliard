/*
 * Landau.cuh
 *
 *  Created on: Nov 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_MODELS_LANDAU_CUH_
#define SRC_CUDA_MODELS_LANDAU_CUH_

#include "../../defs_CUDA.h"

namespace ch {
    void landau_der_bulk_free_energy(double *psi, float *psi_der, int grid_size, float epsilon);
}

#endif /* SRC_CUDA_MODELS_LANDAU_CUH_ */