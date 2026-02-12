/*
 * CUDAIntegrator.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "CUDAIntegrator.h"

namespace ch {

template<int dims>
CUDAIntegrator<dims>::CUDAIntegrator(SimulationState<dims> &sim_state,FreeEnergyModel *model, toml::table &config) : 
        Integrator<dims>(sim_state, model, config) {
    _grid_size = this->_N_bins * model->N_species();

    this->_d_vec_size = this->_N_bins * model->N_species() * sizeof(field_type);
	int d_der_vec_size = this->_N_bins * model->N_species() * sizeof(float);

	this->info("Size of the CUDA direct-space vectors: {} ({} bytes)", this->_N_bins * model->N_species(), this->_d_vec_size);

	this->_h_rho = MultiField<field_type>(this->_N_bins, model->N_species());
	CUDA_SAFE_CALL(cudaMalloc((void **) &sim_state.CUDA_rho, this->_d_vec_size));
    this->_d_rho = sim_state.CUDA_rho;

	CUDA_SAFE_CALL(cudaMalloc((void **) &this->_d_rho_der, d_der_vec_size)); // always float

    this->_h_mobility = MultiField<field_type>(this->_N_bins, model->N_species());
    CUDA_SAFE_CALL(cudaMalloc((void **) &sim_state.CUDA_mobility, this->_d_vec_size));
    this->_d_mobility = sim_state.CUDA_mobility;

    _CPU_GPU();
}

template<int dims>
CUDAIntegrator<dims>::~CUDAIntegrator() {
    CUDA_SAFE_CALL(cudaFree(this->_sim_state.CUDA_rho));
    CUDA_SAFE_CALL(cudaFree(this->_d_rho_der));
    CUDA_SAFE_CALL(cudaFree(this->_sim_state.CUDA_mobility));
}

template<int dims>
void CUDAIntegrator<dims>::_CPU_GPU() {
	for(unsigned int idx = 0; idx < this->_rho.bins(); idx++) {
		for(int species = 0; species < this->_model->N_species(); species++) {
			_h_rho(idx, species) = this->_rho(idx, species);
            _h_mobility(idx, species) = this->_sim_state.mobility(idx, species);
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_rho, _h_rho.data(), _d_vec_size, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(_d_mobility, _h_mobility.data(), _d_vec_size, cudaMemcpyHostToDevice));
}

template<int dims>
void CUDAIntegrator<dims>::sync() {
    if(!_output_ready) {
        CUDA_SAFE_CALL(cudaMemcpy(_h_rho.data(), _d_rho, _d_vec_size, cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaMemcpy(_h_mobility.data(), _d_mobility, _d_vec_size, cudaMemcpyDeviceToHost));

        for(unsigned int idx = 0; idx < this->_rho.bins(); idx++) {
            for(int species = 0; species < this->_model->N_species(); species++) {
                this->_rho(idx, species) = _h_rho(idx, species);
                this->_sim_state.mobility(idx, species) = _h_mobility(idx, species);
            }
        }

        _output_ready = true;
    }
}

template class CUDAIntegrator<1>;
template class CUDAIntegrator<2>;
template class CUDAIntegrator<3>;

} /* namespace ch */
