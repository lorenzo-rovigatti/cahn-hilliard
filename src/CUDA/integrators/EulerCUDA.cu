/*
 * EulerCUDA.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "EulerCUDA.h"

#include "../grid_utils.cuh"

__constant__ int c_N[1]; // number of bins along each direction
__constant__ int c_size[1]; // size of the grid of a single species (N**d)
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1]; // size of the arrays, size * N_species
__constant__ int c_bits[1];

template<int dims> 
__global__ void _Euler_add_surface_term_kernel(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;
    
    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(c_size[0], c_N[0], c_bits[0], rho, IND, dx);
}

template<int dims> 
__global__ void _Euler_integrate_kernel(field_type *rho, float *rho_der, float dx, float dt, float M) {
    if(IND >= c_grid_size[0]) return;

    rho[IND] += (field_type) M * (field_type) _cell_laplacian<dims>(c_size[0], c_N[0], c_bits[0], rho_der, IND, dx) * (field_type) dt;
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
