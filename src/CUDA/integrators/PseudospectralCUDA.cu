/*
 * PseudospectralCUDA.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "PseudospectralCUDA.h"

__constant__ int c_N[1]; // number of bins along each direction
__constant__ int c_size[1]; // size of the grid of a single species (N**d)
__constant__ int c_N_species[1];
__constant__ int c_grid_size[1]; // size of the arrays, size * N_species
__constant__ int c_bits[1];

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

__global__ void _integrate_fft_kernel(cufftFieldComplex *rho_hat, cufftFieldComplex *rho_hat_for_inverse_transform, cufftComplex *f_der_hat, float *sqr_wave_vectors, float *dealiaser, float dt, float M, float k_laplacian, int hat_vector_size) {
    if(IND >= hat_vector_size) return;

    float k2 = sqr_wave_vectors[IND];
    cufftComplex f_der_hat_dealiased = f_der_hat[IND];// * dealiaser[IND];
    
	cufftFieldComplex new_rho_hat = (rho_hat[IND] - dt * M * k2 * f_der_hat_dealiased) / (1.f + dt * M * 2.f * k_laplacian * k2 * k2);
    rho_hat[IND] = new_rho_hat;
    rho_hat_for_inverse_transform[IND] = new_rho_hat / c_size[0];
}

namespace ch {



template<int dims>
PseudospectralCUDA<dims>::PseudospectralCUDA(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
        CUDAIntegrator<dims>(sim_state, model, config) {
	this->_h_rho = MultiField<field_type>(this->_N_bins, model->N_species());
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N, &this->_N_per_dim, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_size, &this->_N_bins, sizeof(int)));
    int N_species = model->N_species();
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_species, &N_species, sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_grid_size, &this->_grid_size, sizeof(int)));
    int bits = (int) std::log2(this->_N_per_dim);
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_bits, &bits, sizeof(int)));

    std::array<int, dims> reciprocal_n;
    reciprocal_n.fill(this->_N_per_dim);
    int hat_grid_size = reciprocal_n[dims - 1] / 2 + 1; 
    for(int i = 0; i < dims - 1; i++) {
        hat_grid_size *= reciprocal_n[i];
    }
    _hat_vector_size = hat_grid_size * model->N_species();

    CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_hat, sizeof(cufftFieldComplex) * _hat_vector_size));
    CUDA_SAFE_CALL(cudaMalloc((void **) &_d_rho_hat_copy, sizeof(cufftFieldComplex) * _hat_vector_size));
    CUDA_SAFE_CALL(cudaMalloc((void **) &_d_f_der_hat, sizeof(cufftComplex) * _hat_vector_size));
    CUDA_SAFE_CALL(cudaMalloc((void **) &_d_sqr_wave_vectors, sizeof(float) * _hat_vector_size));
    CUDA_SAFE_CALL(cudaMalloc((void **) &_d_dealiaser, sizeof(float) * _hat_vector_size));

    std::vector<double> sqr_wave_vectors(_hat_vector_size);
    std::vector<double> dealiaser(_hat_vector_size);

    double nyquist_mode = this->_N_per_dim * M_PI / (this->_N_per_dim * this->_dx) * 2.0 / 3.0;
    if(dims == 1) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(unsigned int i = 0; i < hat_grid_size; i++) {
                double k = 2.0 * M_PI * i / (this->_N_per_dim * this->_dx);
                dealiaser[k_idx] = (k < nyquist_mode) ? 1.0 : 0.0;
                sqr_wave_vectors[k_idx] = SQR(k);
                k_idx++;
            }
        }
    }
    else if(dims == 2) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(int kx_idx = 0; kx_idx < reciprocal_n[0]; kx_idx++) {
                int kx = (kx_idx < reciprocal_n[0] / 2) ? kx_idx : reciprocal_n[0] - kx_idx;
                for(int ky = 0; ky < (reciprocal_n[1] / 2 + 1); ky++) {
                    double k = 2.0 * M_PI * std::sqrt(SQR(kx) + SQR(ky)) / (this->_N_per_dim * this->_dx);
                    dealiaser[k_idx] = (k < nyquist_mode) ? 1.0 : 0.0;
                    sqr_wave_vectors[k_idx] = SQR(k);
                    k_idx++;
                }
            }
        }
    }
    else {
        this->critical("Unsupported number of dimensions {}", dims);
    }

    // sqr_wave_vectors e deliaser are std::vector<double>, so we first have to convert them to std::vector<float>
    // and then we can copy their content to the GPU
    std::vector<float> f_sqr_wave_vectors(sqr_wave_vectors.begin(), sqr_wave_vectors.end());
    CUDA_SAFE_CALL(cudaMemcpy(_d_sqr_wave_vectors, f_sqr_wave_vectors.data(), sizeof(float) * f_sqr_wave_vectors.size(), cudaMemcpyHostToDevice));
    std::vector<float> f_dealiaser(dealiaser.begin(), dealiaser.end());
    CUDA_SAFE_CALL(cudaMemcpy(_d_dealiaser, f_dealiaser.data(), sizeof(float) * f_dealiaser.size(), cudaMemcpyHostToDevice));

    CUFFT_CALL(cufftCreate(&_d_rho_inverse_plan));
    CUFFT_CALL(cufftCreate(&_d_f_der_plan));

#ifdef CUDA_FIELD_FLOAT
    CUFFT_CALL(cufftPlanMany(&_d_rho_inverse_plan, dims, _reciprocal_n.data(), nullptr, 1, odist, nullptr, 1, idist, CUFFT_C2R, model->N_species()));
#else
    CUFFT_CALL(cufftPlanMany(&_d_rho_inverse_plan, dims, reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_Z2D, model->N_species()));
#endif 
    CUFFT_CALL(cufftPlanMany(&_d_f_der_plan, dims, reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_R2C, model->N_species()));

    // Prepare the plan for the forward transform of rho and perform it once
    cufftHandle d_rho_plan;
#ifdef CUDA_FIELD_FLOAT
    CUFFT_CALL(cufftPlanMany(&d_rho_plan, dims, reciprocal_n.data(), nullptr, 1, odist, nullptr, 1, idist, CUFFT_R2C, this->_model->N_species()));
    CUFFT_CALL(cufftExecR2C(d_rho_plan, this->_d_rho, _d_rho_hat));
#else
    CUFFT_CALL(cufftPlanMany(&d_rho_plan, dims, reciprocal_n.data(), nullptr, 1, 0, nullptr, 1, 0, CUFFT_D2Z, this->_model->N_species()));
    CUFFT_CALL(cufftExecD2Z(d_rho_plan, this->_d_rho, _d_rho_hat));
#endif
}

template<int dims>
PseudospectralCUDA<dims>::~PseudospectralCUDA() {
    if(_d_rho_hat != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_rho_hat));
		CUDA_SAFE_CALL(cudaFree(_d_rho_hat_copy));
		CUFFT_CALL(cufftDestroy(_d_rho_inverse_plan));
	}
	if(_d_f_der_hat != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_f_der_hat));
		CUFFT_CALL(cufftDestroy(_d_f_der_plan));
	}
	if(_d_sqr_wave_vectors != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_sqr_wave_vectors));
	}
	if(_d_dealiaser != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_dealiaser));
	}
}

template<int dims>
void PseudospectralCUDA<dims>::evolve() {
    this->_model->der_bulk_free_energy(this->_d_rho, this->_d_rho_der, this->_grid_size);

    CUFFT_CALL(cufftExecR2C(_d_f_der_plan, this->_d_rho_der, _d_f_der_hat));

    const int blocks = this->_grid_size / BLOCK_SIZE + 1;
    _integrate_fft_kernel<<<blocks, BLOCK_SIZE>>>(_d_rho_hat, _d_rho_hat_copy, _d_f_der_hat, 
                                        _d_sqr_wave_vectors, this->_d_dealiaser, this->_dt, this->_M, 
                                        this->_k_laplacian, this->_hat_vector_size);

#ifdef CUDA_FIELD_FLOAT
    CUFFT_CALL(cufftExecC2R(_d_rho_inverse_plan, _d_rho_hat_copy, this->_d_rho));
#else
    CUFFT_CALL(cufftExecZ2D(_d_rho_inverse_plan, _d_rho_hat_copy, this->_d_rho));
#endif

    this->_output_ready = false;
}

template class PseudospectralCUDA<1>;
template class PseudospectralCUDA<2>;

} /* namespace ch */
