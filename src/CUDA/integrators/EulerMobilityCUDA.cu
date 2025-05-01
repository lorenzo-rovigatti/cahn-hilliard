/*
 * EulerMobilityCUDA.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "EulerMobilityCUDA.h"

#include "../grid_utils.cuh"

__constant__ int c_N[1]; // number of bins along each direction
__constant__ int c_size[1]; // size of the grid of a single species (N**d)
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1]; // size of the arrays, size * N_species
__constant__ int c_bits[1];

template<int dims> 
__global__ void _EulerMobility_add_surface_term_kernel(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;
    
    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(c_size[0], c_N[0], c_bits[0], rho, IND, dx);
}

template<int dims>
__global__ void _EulerMobility_compute_flux_kernel(ch::CUDAGrid<dims, ch::CUDAVector<dims>> *grid, field_type *rho, float *rho_der, float dx, float M, float rho_min) {
    if(IND >= grid->total_size) return;

    int species = (IND / grid->species_size);
    int rel_idx = IND % grid->species_size;

    std::array<int, dims> indices;
    for(int d = 0; d < dims; d++) {
        indices[d] = rel_idx % grid->sizes[d];
        rel_idx /= grid->sizes[d];
    }

    ch::CUDAVector<dims> grad;
    for(int d = 0; d < dims; d++) {
        std::array<int, dims> forward = indices;
        forward[d] = (indices[d] + 1) % grid->sizes[d];
        int idx_p = grid->index(forward, species);

        field_type M_idx = M * rho[IND] / (rho[IND] + rho_min);
        field_type M_p = M * rho[idx_p] / (rho[idx_p] + rho_min);
        field_type M_flux = 0.5f * (M_idx + M_p);

        grad[d] = M_flux * (rho_der[idx_p] - rho_der[IND]) / dx;
    }

    grid->data()[IND] = grad;
}

template<int dims> 
__global__ void _EulerMobility_integrate_kernel(field_type *rho, ch::CUDAGrid<dims, ch::CUDAVector<dims>> *grid, float dx, float dt) {
    if(IND >= grid->total_size) return;

    int species = (IND / grid->species_size);
    int rel_idx = IND % grid->species_size;

    std::array<int, dims> indices;
    for(int d = 0; d < dims; d++) {
        indices[d] = rel_idx % grid->sizes[d];
        rel_idx /= grid->sizes[d];
    }

    field_type divergence = 0.f;
    for(int d = 0; d < dims; d++) {
        std::array<int, dims> backward = indices;
        backward[d] = (indices[d] - 1 + grid->sizes[d]) % grid->sizes[d];

        int idx_m = grid->index(backward, species);
        divergence += (grid->data()[IND][d] - grid->data()[idx_m][d]) / dx;
    }

    rho[IND] += divergence * (field_type) dt;
}

namespace ch {

template<int dims>
EulerMobilityCUDA<dims>::EulerMobilityCUDA(FreeEnergyModel *model, toml::table &config) : CUDAIntegrator<dims>(model, config) {
    this->_d_vec_size = this->_N_bins * model->N_species() * sizeof(field_type);
	int d_der_vec_size = this->_N_bins * model->N_species() * sizeof(float);

	this->info("Size of the CUDA direct-space vectors: {} ({} bytes)", this->_N_bins * model->N_species(), this->_d_vec_size);

    _rho_min = this->template _config_value<float>(config, "mobility.rho_min");

    int N_species = model->N_species();
    CUDAGrid<dims, CUDAVector<dims>> h_flux(this->_N_per_dim, N_species);
    CUDA_SAFE_CALL(cudaMalloc((CUDAGrid<dims, CUDAVector<dims>> **) &_d_flux, sizeof(CUDAGrid<dims, CUDAVector<dims>>)));
    CUDA_SAFE_CALL(cudaMemcpy(_d_flux, &h_flux, sizeof(CUDAGrid<dims, CUDAVector<dims>>), cudaMemcpyHostToDevice));

	this->_h_rho = RhoMatrix<field_type>(this->_N_bins, N_species);
	CUDA_SAFE_CALL(cudaMalloc((void **) &this->_d_rho, this->_d_vec_size));
	CUDA_SAFE_CALL(cudaMalloc((void **) &this->_d_rho_der, d_der_vec_size)); // always float

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &this->_N_per_dim, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &this->_N_bins, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &this->_grid_size, sizeof(int)));
    int bits = (int) std::log2(this->_N_per_dim);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_bits, &bits, sizeof(int)));
}

template<int dims>
EulerMobilityCUDA<dims>::~EulerMobilityCUDA() {
    if(this->_d_rho != nullptr) {
		CUDA_SAFE_CALL(cudaFree(this->_d_rho));
	}
	if(this->_d_rho_der != nullptr) {
		CUDA_SAFE_CALL(cudaFree(this->_d_rho_der));
	}
    if(_d_flux != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_flux));
    }
}

template<int dims>
void EulerMobilityCUDA<dims>::evolve() {
    this->_model->der_bulk_free_energy(this->_d_rho, this->_d_rho_der, this->_grid_size);

    const int blocks = this->_grid_size / BLOCK_SIZE + 1;
    _EulerMobility_add_surface_term_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, this->_d_rho_der, this->_dx, this->_k_laplacian);
    _EulerMobility_compute_flux_kernel<dims><<<blocks, BLOCK_SIZE>>>(_d_flux, this->_d_rho, this->_d_rho_der, this->_dx, this->_M, _rho_min);
    _EulerMobility_integrate_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, _d_flux, this->_dx, this->_dt);

    this->_output_ready = false;
}

template class EulerMobilityCUDA<1>;
template class EulerMobilityCUDA<2>;

} /* namespace ch */
