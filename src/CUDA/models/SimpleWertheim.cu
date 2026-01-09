#include "SimpleWertheim.cuh"

#include <stdio.h>

__device__ float _CUDA_X(float rho, float two_valence_delta) {
	float sqrt_argument = 2.f * two_valence_delta * rho;
	return (sqrt_argument < 1e-3f) ? 1.f - sqrt_argument / 2.f : (-1.f + sqrtf(1.f + 2.f * two_valence_delta * rho)) / (two_valence_delta * rho);
}

__global__ void _compute_mobility(field_type *rho, field_type *mobility, float M0, int grid_size, int valence, float two_valence_delta) {
	if(IND >= grid_size) return;

	float rho_tot = rho[IND];
	float X = _CUDA_X(rho_tot, two_valence_delta);
	float D = M0 * powf(X, (float) valence);

	mobility[IND] = rho_tot * D;
}

__global__ void _compute_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2, int valence, float two_valence_delta) {
	if(IND >= grid_size) return;

	float rho_ind = rho[IND];
	float der_f_ref = logf(rho_ind) + 2.f * B2 * rho_ind;
	float X = _CUDA_X(rho_ind, two_valence_delta);
	float der_f_bond = valence * logf(X);

	rho_der[IND] = der_f_ref + der_f_bond;
}

namespace ch {
	void simple_wertheim_mobility(field_type *rho, field_type *mobility, double M0, int grid_size, int valence, float two_valence_delta) {
		const int blocks = grid_size / BLOCK_SIZE + 1;
		_compute_mobility<<<blocks, BLOCK_SIZE>>>(rho, mobility, M0, grid_size, valence, two_valence_delta);
	}

    void simple_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size, float B2, int valence, float two_valence_delta) {
        const int blocks = grid_size / BLOCK_SIZE + 1;
        _compute_der_bulk_free_energy<<<blocks, BLOCK_SIZE>>>(rho, rho_der, grid_size, B2, valence, two_valence_delta);
    }

}
