#include "CahnHilliard.cuh"

__constant__ int c_N[1];
__constant__ int c_size[1];
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1];

template<int dims> 
__device__ float _cell_laplacian(number *rho, int idx, float dx);


template<> 
__device__ float _cell_laplacian<1>(number *rho, int idx, float dx) {
    int N_minus_one = c_N[0] - 1;
    int idx_m = (idx - 1 + c_N[0]) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (rho[idx_m] + rho[idx_p] - 2.f * rho[idx]) / (dx * dx);
}

template<> 
__device__ float _cell_laplacian<2>(number *rho, int idx, float dx) {
    int N_minus_one = c_N[0];
    int idx_m = (idx - 1 + c_N[0]) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (rho[idx_m] + rho[idx_p] - 2.f * rho[idx]) / (dx * dx);
}

template<int dims> 
__global__ void _add_surface_term(number *rho, number *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;

    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(rho, IND, dx);
}

template<int dims> 
__global__ void _integrate(number *rho, number *rho_der, float dx, float dt, float M) {
    if(IND >= c_grid_size[0]) return;

    rho[IND] += M * _cell_laplacian<dims>(rho_der, IND, dx) * dt;
}

namespace ch {

int grid_size;
void init_symbols(int N, int size, int N_species) {
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &N, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &size, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    grid_size = size * N_species;
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &grid_size, sizeof(int)));
}

template<int dims> 
void add_surface_term<dims>(number *rho, number *rho_der, number dx, number k_laplacian) {
    const int blocks = grid_size / BLOCK_SIZE + 1;
    _add_surface_term<dims><<<blocks, BLOCK_SIZE>>>(rho, rho_der, dx, k_laplacian);
}

template<int dims> 
void integrate<dims>(number *rho, number *rho_der, number dx, number dt, number M) {
    const int blocks = grid_size / BLOCK_SIZE + 1;
    _integrate<dims><<<blocks, BLOCK_SIZE>>>(rho, rho_der, dx, dt, M);
}

template void add_surface_term<1>(number *rho, number *rho_der, number dx, number k_laplacian);
template void add_surface_term<2>(number *rho, number *rho_der, number dx, number k_laplacian);

template void integrate<1>(number *rho, number *rho_der, number dx, number dt, number M);
template void integrate<2>(number *rho, number *rho_der, number dx, number dt, number M);

}
