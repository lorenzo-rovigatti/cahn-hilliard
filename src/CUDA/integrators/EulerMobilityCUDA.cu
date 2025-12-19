/*
 * EulerMobilityCUDA.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "EulerMobilityCUDA.h"

#include "../grid_utils.cuh"
#include "../../utils/utility_functions.h"

__constant__ int c_N[1]; // number of bins along each direction
__constant__ int c_size[1]; // size of the grid of a single species (N**d)
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1]; // size of the arrays, size * N_species
__constant__ int c_bits[1];
__constant__ float c_noise_factor[1];

template<int dims> 
__global__ void _EulerMobility_add_surface_term_kernel(field_type *rho, float *rho_der, float dx, float k_laplacian) {
    if(IND >= c_grid_size[0]) return;
    
    rho_der[IND] -= 2.f * k_laplacian * _cell_laplacian<dims>(c_size[0], c_N[0], c_bits[0], rho, IND, dx);
}

template<int dims>
__global__ void _EulerMobility_compute_flux_kernel(
    ch::CUDAGrid<dims, ch::CUDAVector<dims>> *grid, 
    field_type *rho, 
    float *rho_der,
    curandState *rand_states,
    float dx, 
    float M, 
    float rho_min,
    bool with_noise
) {
    if(IND >= grid->total_size) return;

    int species = (IND / grid->species_size);
    int rel_idx = IND % grid->species_size;

    std::array<int, dims> indices;
    for(int d = 0; d < dims; d++) {
        indices[d] = rel_idx % grid->sizes[d];
        rel_idx /= grid->sizes[d];
    }

    ch::CUDAVector<dims> grad;

    curandState local_state;
    if(with_noise) {
        local_state = rand_states[IND];
    }

    for(int d = 0; d < dims; d++) {
        std::array<int, dims> forward = indices;
        forward[d] = (indices[d] + 1) % grid->sizes[d];
        int idx_p = grid->index(forward, species);

        field_type M_idx = M * rho[IND] / (rho[IND] + rho_min);
        field_type M_p = M * rho[idx_p] / (rho[idx_p] + rho_min);
        field_type M_flux = 0.5f * (M_idx + M_p);

        field_type deterministic = M_flux * (rho_der[idx_p] - rho_der[IND]) / dx;
        field_type stochastic = 0.0f;
        if(with_noise) {
            field_type noise_amplitude = sqrtf(M_flux) * c_noise_factor[0];
            field_type gaussian = curand_normal(&local_state);
            stochastic = noise_amplitude * gaussian;
        }

        grad[d] = deterministic + stochastic;
    }

    if(with_noise) {
        rand_states[IND] = local_state;
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

__global__ void init_rand_states(curandState *states, unsigned long seed, int total_size) {
    if(IND >= total_size) return;

    // Each thread gets a unique sequence
    curand_init(seed, IND, 0, &states[IND]);
}

namespace ch {

template<int dims>
EulerMobilityCUDA<dims>::EulerMobilityCUDA(SimulationState &sim_state,FreeEnergyModel *model, toml::table &config) : 
        CUDAIntegrator<dims>(sim_state, model, config) {
    _rho_min = this->template _config_value<float>(config, "mobility.rho_min");

    int N_species = model->N_species();
    _h_flux = new CUDAGrid<dims, CUDAVector<dims>>(this->_N_per_dim, N_species);
    CUDA_SAFE_CALL(cudaMalloc((CUDAGrid<dims, CUDAVector<dims>> **) &_d_flux, sizeof(CUDAGrid<dims, CUDAVector<dims>>)));
    CUDA_SAFE_CALL(cudaMemcpy(_d_flux, _h_flux, sizeof(CUDAGrid<dims, CUDAVector<dims>>), cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &this->_N_per_dim, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &this->_N_bins, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &this->_grid_size, sizeof(int)));
    int bits = (int) std::log2(this->_N_per_dim);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_bits, &bits, sizeof(int)));

    _with_noise = this->template _config_optional_value<bool>(config, "mobility.with_noise", false);
	if(_with_noise) {
		double noise_rescale_factor = this->template _config_optional_value<double>(config, "mobility.noise_rescale_factor", 1.0);
		double noise_factor = std::sqrt(2.0 / (pow_dims<dims>(this->_dx) * this->_dt)) * noise_rescale_factor;
        COPY_NUMBER_TO_FLOAT(c_noise_factor, noise_factor);

        CUDA_SAFE_CALL(cudaMalloc(&_d_rand_states, _h_flux->total_size * sizeof(curandState)));
		long long int seed = this->template _config_optional_value<long long int>(config, "seed", std::time(NULL));
        int gridSize = (_h_flux->total_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
        init_rand_states<<<gridSize, BLOCK_SIZE>>>(_d_rand_states, seed, _h_flux->total_size);
        cudaDeviceSynchronize();

        this->info("Integrating the Cahn-Hilliard equation on CUDA with non-constant mobility and noise (noise_factor = {})", noise_factor);
	}
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
    if(_d_rand_states != nullptr) {
        CUDA_SAFE_CALL(cudaFree(_d_rand_states));
    }
}

template<int dims>
void EulerMobilityCUDA<dims>::evolve() {
    this->_model->der_bulk_free_energy(this->_d_rho, this->_d_rho_der, this->_grid_size);

    const int blocks = this->_grid_size / BLOCK_SIZE + 1;
    _EulerMobility_add_surface_term_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, this->_d_rho_der, this->_dx, this->_k_laplacian);
    _EulerMobility_compute_flux_kernel<dims><<<blocks, BLOCK_SIZE>>>(_d_flux, this->_d_rho, this->_d_rho_der, this->_d_rand_states, this->_dx, this->_M, _rho_min, _with_noise);
    _EulerMobility_integrate_kernel<dims><<<blocks, BLOCK_SIZE>>>(this->_d_rho, _d_flux, this->_dx, this->_dt);

    this->_output_ready = false;
}

template class EulerMobilityCUDA<1>;
template class EulerMobilityCUDA<2>;

} /* namespace ch */
