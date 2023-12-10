#include "CahnHilliard.cuh"

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
__global__ void _add_surface_term(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;
    
    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(rho, IND, dx);
}

template<int dims> 
__global__ void _integrate(field_type *rho, float *rho_der, float dx, float dt, float M) {
    if(IND >= c_grid_size[0]) return;

    rho[IND] += (field_type) M * (field_type) _cell_laplacian<dims>(rho_der, IND, dx) * (field_type) dt;
}

__device__ cufftComplex operator*(const cufftComplex &lhs, const float &rhs) {
    return cufftComplex({lhs.x * rhs, lhs.y * rhs});
}

__device__ cufftComplex operator*(const float &lhs, const cufftComplex &rhs) {
    return rhs * lhs;
}

__device__ cufftFieldComplex operator-(const cufftFieldComplex &lhs, const cufftComplex &rhs) {
    return cufftFieldComplex({lhs.x - rhs.x, lhs.y - rhs.y});
}

__device__ cufftFieldComplex operator/(const cufftFieldComplex &lhs, const float &rhs) {
    return cufftFieldComplex({lhs.x / rhs, lhs.y / rhs});
}

__global__ void _integrate_fft(cufftFieldComplex *rho_hat, cufftFieldComplex *rho_hat_for_inverse_transform, cufftComplex *f_der_hat, float *sqr_wave_vectors, float *dealiaser, float dt, float M, float k_laplacian, int size_hat) {
    if(IND >= size_hat) return;

    float k2 = sqr_wave_vectors[IND];
    cufftComplex f_der_hat_dealiased = f_der_hat[IND];// * dealiaser[IND];
    
	cufftFieldComplex new_rho_hat = (rho_hat[IND] - dt * M * k2 * f_der_hat_dealiased) / (1.f + dt * M * 2.f * k_laplacian * k2 * k2);
    rho_hat[IND] = new_rho_hat;
    rho_hat_for_inverse_transform[IND] = new_rho_hat / c_size[0];
}

namespace ch {

int grid_size;
void init_symbols(int N, int size, int N_species) {
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &N, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &size, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    grid_size = size * N_species;
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &grid_size, sizeof(int)));
    int bits = (int) std::log2(N);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_bits, &bits, sizeof(int)));
}

template<int dims> 
void add_surface_term<dims>(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    const int blocks = grid_size / BLOCK_SIZE + 1;
    _add_surface_term<dims><<<blocks, BLOCK_SIZE>>>(rho, rho_der, dx, k_laplacian);
}

template<int dims> 
void integrate<dims>(field_type *rho, float *rho_der, float dx, float dt, float M) {
    const int blocks = grid_size / BLOCK_SIZE + 1;
    _integrate<dims><<<blocks, BLOCK_SIZE>>>(rho, rho_der, dx, dt, M);
}

void integrate_fft(cufftFieldComplex *rho_hat, cufftFieldComplex *rho_hat_for_inverse_transform, cufftComplex *f_der_hat, float *sqr_wave_vectors, float *dealiaser, float dt, float M, float k_laplacian, int hat_size) {
    const int blocks = hat_size / BLOCK_SIZE + 1;
    _integrate_fft<<<blocks, BLOCK_SIZE>>>(rho_hat, rho_hat_for_inverse_transform, f_der_hat, sqr_wave_vectors, dealiaser, dt, M, k_laplacian, hat_size);
}

template void add_surface_term<1>(field_type *rho, float *rho_der, float dx, float k_laplacian);
template void add_surface_term<2>(field_type *rho, float *rho_der, float dx, float k_laplacian);

template void integrate<1>(field_type *rho, float *rho_der, float dx, float dt, float M);
template void integrate<2>(field_type *rho, float *rho_der, float dx, float dt, float M);

}
