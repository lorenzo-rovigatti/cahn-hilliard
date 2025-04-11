/*
 * GenericWertheim.cuh
 *
 *  Created on: Apr 11, 2025
 *      Author: lorenzo
 */

#ifndef SRC_CUDA_MODELS_GENERICWERTHEIM_CUH_
#define SRC_CUDA_MODELS_GENERICWERTHEIM_CUH_


#include <vector>

#include "../../defs_CUDA.h"

#define MAX_PATCH_INTERACTIONS 10

namespace ch {
    void init_generic_wertheim_symbols(
        int N_species, int N_patches, std::vector<int> &unique_patch_ids, std::vector<ushort2> &species_patches, 
        std::vector<int> &interacting_species, std::vector<ushort2> &interacting_patches, std::vector<double> &B2, 
        std::vector<double> &delta);
    void generic_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size);
}

#endif /* SRC_CUDA_MODELS_GENERICWERTHEIM_CUH_ */