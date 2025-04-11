#include "GenericWertheim.cuh"

#include <stdio.h>

__constant__ int c_N_species[1];
__constant__ int c_N_patches[1];
__constant__ int c_N_unique_patches[1];

__global__ void _generic_wertheim_update_X(field_type *rho, int grid_size, int *unique_patch_ids, 
	int *interacting_species, ushort2 *interacting_patches, float *delta, float *X) {

	if(IND >= grid_size) return;

	float tolerance = 1e-6f;
	int max_iter = 10000;

	for(int iter = 0; iter < max_iter; iter++) {
		float max_delta = 0.f;

		for(int i = 0; i < c_N_unique_patches[0]; i++) {
			int patch = unique_patch_ids[i];
			float sum = 0.f;
			for(int j = 0; j < MAX_PATCH_INTERACTIONS; j++) {
				int species = interacting_species[patch * MAX_PATCH_INTERACTIONS + j];
				if(species != -1) {
					float rho_species = rho[species * grid_size + IND];
					for(int k = 0; k < c_N_patches[0]; k++) {
						ushort2 other_patch = interacting_patches[patch * c_N_patches[0] * MAX_PATCH_INTERACTIONS + j * c_N_patches[0] + k];
						float int_delta = delta[patch * c_N_patches[0] + other_patch.x];
						// printf("%d %d %d %d %lf\n", patch, species, other_patch.x, other_patch.y, int_delta);
						float patch_X = X[other_patch.x * grid_size + IND];
						sum += other_patch.y * rho_species * patch_X * int_delta;
					}
				}
			}

			float new_X = 1.f / (1.f + sum);
			float patch_X = X[patch * grid_size + IND];
			max_delta = max(max_delta, fabs(new_X - patch_X));
			X[patch * grid_size + IND] = new_X;
		}

		if(max_delta < tolerance) {
			return;
		}
	}
}

__global__ void _compute_generic_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size, float *B2, ushort2 *species_patches, float *X) {
	if(IND >= vec_size) return;

	int grid_size = vec_size / c_N_species[0];
	int species = IND / grid_size;
    int rel_idx = IND % grid_size;

	float rho_species = (float) rho[species * grid_size + rel_idx];
	float B2_contrib = 0.f;
	for(int other_species = 0; other_species < c_N_species[0]; other_species++) {
		float rho_other_species = (float) rho[other_species * grid_size + rel_idx];
		B2_contrib += 2.f * rho_other_species * B2[species * c_N_species[0] + other_species];
	}
	float der_f_ref = logf(rho_species) + B2_contrib;

	float der_f_bond = 0.f;
	for(int p = 0; p < c_N_patches[0]; p++) {
		ushort2 patch = species_patches[species * c_N_patches[0] + p];
		float patch_X = X[patch.x * grid_size + rel_idx];
		der_f_bond += patch.y * logf(patch_X); 
	}

	rho_der[IND] = der_f_ref + der_f_bond;
}

namespace ch {
	int *d_unique_patch_ids;
	float *d_B2;
	float *d_delta;
	ushort2 *d_species_patches;
	int *d_interacting_species;
	ushort2 *d_interacting_patches;
	float *d_X = nullptr;
	int h_N_species, h_N_patches; // needed to allocate d_X

	void init_generic_wertheim_symbols(
		int N_species, int N_patches, std::vector<int> &unique_patch_ids, std::vector<ushort2> &species_patches, 
		std::vector<int> &interacting_species, std::vector<ushort2> &interacting_patches, std::vector<double> &B2, 
		std::vector<double> &delta) {

		h_N_species = N_species;
		h_N_patches = N_patches;
		int N_unique_patches = unique_patch_ids.size();

		CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_patches, &N_patches, sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_unique_patches, &N_unique_patches, sizeof(int)));

		CUDA_SAFE_CALL(cudaMalloc((void **) &d_unique_patch_ids, unique_patch_ids.size() * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpy(d_unique_patch_ids, unique_patch_ids.data(), unique_patch_ids.size() * sizeof(int), cudaMemcpyHostToDevice));

		std::vector<float> B2_float(B2.begin(), B2.end());
		CUDA_SAFE_CALL(cudaMalloc((void **) &d_B2, B2_float.size() * sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpy(d_B2, B2_float.data(), B2_float.size() * sizeof(float), cudaMemcpyHostToDevice));

		std::vector<float> delta_float(delta.begin(), delta.end());
		CUDA_SAFE_CALL(cudaMalloc((void **) &d_delta, delta_float.size() * sizeof(float)));
		CUDA_SAFE_CALL(cudaMemcpy(d_delta, delta_float.data(), delta_float.size() * sizeof(float), cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(cudaMalloc((void **) &d_species_patches, species_patches.size() * sizeof(ushort2)));
		CUDA_SAFE_CALL(cudaMemcpy(d_species_patches, species_patches.data(), species_patches.size() * sizeof(ushort2), cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(cudaMalloc((void **) &d_interacting_species, interacting_species.size() * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpy(d_interacting_species, interacting_species.data(), interacting_species.size() * sizeof(int), cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(cudaMalloc((void **) &d_interacting_patches, interacting_patches.size() * sizeof(ushort2)));
		CUDA_SAFE_CALL(cudaMemcpy(d_interacting_patches, interacting_patches.data(), interacting_patches.size() * sizeof(ushort2), cudaMemcpyHostToDevice));
	}

    void generic_wertheim_der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
		int grid_size = vec_size / h_N_species;
		// we allocate d_X the first time this function is called
		if(d_X == nullptr) {
			int X_size = grid_size * h_N_patches;
			CUDA_SAFE_CALL(cudaMalloc((void **) &d_X, X_size * sizeof(float)));
			CUDA_SAFE_CALL(cudaMemset(d_X, 0, X_size));
		}

		// one thread for each spatial bin of the X grid, which has dimensions "bins * N_patches
		const int X_blocks = grid_size / BLOCK_SIZE + 1;
		_generic_wertheim_update_X<<<X_blocks, BLOCK_SIZE>>>(rho, grid_size, d_unique_patch_ids, d_interacting_species, d_interacting_patches, d_delta, d_X);

		// one thread for each element of the rho grid, which has dimensions "bins * N_species"
        const int c_blocks = vec_size / BLOCK_SIZE + 1;
        _compute_generic_wertheim_der_bulk_free_energy<<<c_blocks, BLOCK_SIZE>>>(rho, rho_der, vec_size, d_B2, d_species_patches, d_X);
    }

}
