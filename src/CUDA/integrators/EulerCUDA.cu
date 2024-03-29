/*
 * EulerCUDA.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "EulerCUDA.h"

__constant__ int c_N[1]; // number of bins along each direction
__constant__ int c_size[1]; // size of the grid of a single species (N**d)
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1]; // size of the arrays, size * N_species
__constant__ int c_bits[1];

template<int dims> 
__device__ void _fill_coords(int coords[dims], int idx) {
    for(int d = 0; d < dims; d++) {
		coords[d] = idx & (c_N[0] - 1);
		idx >>= c_bits[0]; // divide by N
	}
}

template<int dims>
__device__ int _cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= c_bits[0]; // multiply by N
	}
	return idx;
}

template<int dims, typename number> 
__device__ float _cell_laplacian(number *field, int idx, float dx) {
    int base_idx = (idx / c_size[0]) * c_size[0];
    int rel_idx = idx % c_size[0];
    int N_minus_one = c_N[0] - 1;

    if(dims == 1) {
        int rel_idx_m = (rel_idx - 1 + c_N[0]) & N_minus_one;
        int rel_idx_p = (rel_idx + 1) & N_minus_one;

        return ((float) field[base_idx + rel_idx_m] + (float) field[base_idx + rel_idx_p] - 2.f * (float) field[idx]) / (dx * dx);
    }
    else if(dims == 2) {
        int coords_xy[2];
        _fill_coords<2>(coords_xy, rel_idx);

        int coords_xmy[2] = {
                (coords_xy[0] - 1 + c_N[0]) & N_minus_one,
                coords_xy[1]
        };

        int coords_xym[2] = {
                coords_xy[0],
                (coords_xy[1] - 1 + c_N[0]) & N_minus_one
        };

        int coords_xpy[2] = {
                (coords_xy[0] + 1) & N_minus_one,
                coords_xy[1]
        };

        int coords_xyp[2] = {
                coords_xy[0],
                (coords_xy[1] + 1) & N_minus_one
        };

        return (
                (float) field[base_idx + _cell_idx<2>(coords_xmy)] +
                (float) field[base_idx + _cell_idx<2>(coords_xpy)] +
                (float) field[base_idx + _cell_idx<2>(coords_xym)] +
                (float) field[base_idx + _cell_idx<2>(coords_xyp)] -
                4.f * (float) field[idx])
                / (dx * dx);
    }
}

template<int dims> 
__global__ void _Euler_add_surface_term_kernel(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;
    
    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(rho, IND, dx);
}

template<int dims> 
__global__ void _Euler_integrate_kernel(field_type *rho, float *rho_der, float dx, float dt, float M) {
    if(IND >= c_grid_size[0]) return;

    rho[IND] += (field_type) M * (field_type) _cell_laplacian<dims>(rho_der, IND, dx) * (field_type) dt;
}

namespace ch {



template<int dims>
EulerCUDA<dims>::EulerCUDA(FreeEnergyModel *model, toml::table &config) : CUDAIntegrator<dims>(model, config) {
    this->_d_vec_size = this->_N_bins * model->N_species() * sizeof(field_type);
	int d_der_vec_size = this->_N_bins * model->N_species() * sizeof(float);

	this->info("Size of the CUDA direct-space vectors: {} ({} bytes)", this->_N_bins * model->N_species(), this->_d_vec_size);

	this->_h_rho = RhoMatrix<field_type>(this->_N_bins, model->N_species());
	CUDA_SAFE_CALL(cudaMalloc((void **) &this->_d_rho, this->_d_vec_size));
	CUDA_SAFE_CALL(cudaMalloc((void **) &this->_d_rho_der, d_der_vec_size)); // always float

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &this->_N_per_dim, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &this->_N_bins, sizeof(int)));
    int N_species = model->N_species();
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &this->_grid_size, sizeof(int)));
    int bits = (int) std::log2(this->_N_per_dim);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_bits, &bits, sizeof(int)));
}

template<int dims>
EulerCUDA<dims>::~EulerCUDA() {
    if(this->_d_rho != nullptr) {
		CUDA_SAFE_CALL(cudaFree(this->_d_rho));
	}
	if(this->_d_rho_der != nullptr) {
		CUDA_SAFE_CALL(cudaFree(this->_d_rho_der));
	}
}

template<int dims>
void EulerCUDA<dims>::evolve() {
    this->_model->der_bulk_free_energy(this->_d_rho, this->_d_rho_der, this->_grid_size);

    const int blocks = this->_grid_size / BLOCK_SIZE + 1;
    _Euler_add_surface_term_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, this->_d_rho_der, this->_dx, this->_k_laplacian);
    _Euler_integrate_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, this->_d_rho_der, this->_dx, this->_dt, this->_M);

    this->_output_ready = false;
}

template class EulerCUDA<1>;
template class EulerCUDA<2>;

} /* namespace ch */
