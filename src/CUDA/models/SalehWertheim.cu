#include "SalehWertheim.cuh"

#include <stdio.h>

#define N_SPECIES 3

__constant__ int c_valence[2];
__constant__ int c_linker_half_valence[1];
__constant__ float c_delta_AA[1];
__constant__ float c_delta_BB[1];

__device__ float _der_contribution(float rhos[3], int species) {
	float delta = (species == 0) ? c_delta_AA[0] : c_delta_BB[0];
	float rho_factor =  delta * (c_valence[species] * rhos[species] + c_linker_half_valence[0] * rhos[2]);
	float X = (-1.f + sqrtf(1.f + 4.f * rho_factor)) / (2.f * rho_factor);
	return logf(X);
}

__global__ void _compute_saleh_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2) {
	if(IND >= grid_size) return;

	int size = grid_size / N_SPECIES;
	int species = IND / size;
    int rel_idx = IND % size;

	float rhos[3] = {
		(float) rho[rel_idx],
		(float) rho[size + rel_idx],
		(float) rho[2 * size + rel_idx]
	};

    // the ideal + B2 part is computed analytically
	float der_f_ref = logf(rhos[species]);
	for(int i = 0; i < N_SPECIES; i++) {
		der_f_ref += 2.f * B2 * rhos[i];
	}

	float der_f_bond;
	if(species  < 2) {
		der_f_bond = c_valence[species] * _der_contribution(rhos, species);
	}
	else {
		der_f_bond = c_linker_half_valence[0] * (_der_contribution(rhos, 0) + _der_contribution(rhos, 1));
	}

	rho_der[IND] = der_f_ref + der_f_bond;
}

namespace ch {
	void init_saleh_symbols(std::vector<int> &valence, int linker_half_valence, float delta_AA, float delta_BB) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_valence, valence.data(), 2 * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_linker_half_valence, &linker_half_valence, sizeof(int)));
		COPY_NUMBER_TO_FLOAT(c_delta_AA, delta_AA);
		COPY_NUMBER_TO_FLOAT(c_delta_BB, delta_BB);
	}

    void saleh_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2) {
        const int blocks = grid_size / BLOCK_SIZE + 1;
        _compute_saleh_der_bulk_free_energy<<<blocks, BLOCK_SIZE>>>(rho, rho_der, grid_size, B2);
    }

}
