#include "RicciWertheim.cuh"

#include <stdio.h>

#define N_SPECIES 2
#define SQR(X) ((X) * (X))

__constant__ float c_delta_00[1];
__constant__ float c_delta_12[1];

__global__ void _compute_ricci_der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size, float B2) {
	if(IND >= vec_size) return;

	int grid_size = vec_size / N_SPECIES;
	int species = IND / grid_size;
    int rel_idx = IND % grid_size;

	float rhos[2] = {
		(float) rho[rel_idx],
		(float) rho[grid_size + rel_idx],
	};

	float rho_tot = rhos[0] + rhos[1];
	float der_f_ref = logf(rhos[species]) + 2 * B2 * rho_tot;

	float der_f_bond;
	float A = rhos[1] * c_delta_12[0];
	float B = rhos[0] * c_delta_12[0];
	if(species == 0) {
		float X_0 = (-1 + std::sqrt((1 + 12 * rhos[0] * c_delta_00[0]))) / (6 * rhos[0] * c_delta_00[0]);
		float X_1 = (B - 1 - 2 * A + std::sqrt(SQR(1 + 2 * A - B) + 4 * B)) / (2 * B);
		der_f_bond = 3 * logf(X_0) + logf(X_1);
	}
	else {
		float X_2 = (-B + 2 * A - 1 + sqrtf(SQR(B - 2 * A + 1) + 8 * A)) / (4 * A);
		der_f_bond = 2 * logf(X_2);
	}

	rho_der[IND] = der_f_ref + der_f_bond;
}

namespace ch {
	void init_ricci_symbols(float delta_00, float delta_12) {
		COPY_NUMBER_TO_FLOAT(c_delta_00, delta_00);
		COPY_NUMBER_TO_FLOAT(c_delta_12, delta_12);
	}

    void ricci_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size, float B2) {
        const int blocks = vec_size / BLOCK_SIZE + 1;
        _compute_ricci_der_bulk_free_energy<<<blocks, BLOCK_SIZE>>>(rho, rho_der, vec_size, B2);
    }

}
