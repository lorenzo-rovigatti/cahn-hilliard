#include "Landau.cuh"

#include <stdio.h>

__global__ void _compute_landau_der_bulk_free_energy(field_type *psi, float *psi_der, int grid_size, float epsilon) {
	if(IND >= grid_size) return;

    float op = psi[IND];
	psi_der[IND] = -epsilon * op + op * op * op;
}

namespace ch {

    void landau_der_bulk_free_energy(field_type *psi, float *psi_der, int grid_size, float epsilon) {
        const int blocks = grid_size / BLOCK_SIZE + 1;
        _compute_landau_der_bulk_free_energy<<<blocks, BLOCK_SIZE>>>(psi, psi_der, grid_size, epsilon);
    }

}
