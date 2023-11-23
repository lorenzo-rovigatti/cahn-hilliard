#include "SalehWertheim.cuh"

#include <stdio.h>

#define N_SPECIES 3

__constant__ int c_valence[2];
__constant__ int c_linker_half_valence[1];
__constant__ float c_delta_AA[1];
__constant__ float c_delta_BB[1];

__device__ double _saleh_bonding_free_energy(double rhos[3]) {
	double rho_factor =  c_delta_AA[0] * (c_valence[0] * rhos[0] + c_linker_half_valence[0] * rhos[2]);
	double X_1A = (-1.0 + sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_1 = log(X_1A) - X_1A / 2.0 + 0.5;

	rho_factor =  c_delta_BB[0] * (c_valence[1] * rhos[1] + c_linker_half_valence[0] * rhos[2]);
	double X_2B = (-1.0 + sqrt(1.0 + 4.0 * rho_factor)) / (2.0 * rho_factor);
	double fe_part_2 = log(X_2B) - X_2B / 2.0 + 0.5;

	double bonding_fe =
			rhos[0] * c_valence[0] * fe_part_1 +
			rhos[1] * c_valence[1] * fe_part_2 +
			rhos[2] * c_linker_half_valence[0] * (fe_part_1 + fe_part_2);

	return bonding_fe;
}

__global__ void _compute_saleh_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2) {
	if(IND >= grid_size) return;

	int size = grid_size / N_SPECIES;
	int species = IND / size;
    int rel_idx = IND % size;

	double rhos[3] = {
		rho[rel_idx],
		rho[size + rel_idx],
		rho[2 * size + rel_idx]
	};

    // the ideal + B2 part is computed analytically
	float der_f_ref = logf(rhos[species]);
	for(int i = 0; i < N_SPECIES; i++) {
		der_f_ref += 2.f * B2 * rhos[i];
	}

	// the bonding part is computed numerically
	double delta_rho_i = rhos[species] * 1e-5;
	double fe_r = _saleh_bonding_free_energy(rhos);
	rhos[species] += delta_rho_i;
	double fe_rdr = _saleh_bonding_free_energy(rhos);
	double der_f_bond = (fe_rdr - fe_r) / delta_rho_i;

	// if(rel_idx == 0 && species == 1) printf("%d %lf\n", IND, der_f_ref + der_f_bond);

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
