/*
 * CUDAIntegrator.cu
 *
 * Created on: 3/29/2024
 *     Author: Lorenzo
*/

#include "CUDAIntegrator.h"

namespace ch {

template<int dims>
CUDAIntegrator<dims>::CUDAIntegrator(FreeEnergyModel *model, toml::table &config) : Integrator<dims>(model, config) {
    _grid_size = this->_N_bins * model->N_species();
}

template<int dims>
CUDAIntegrator<dims>::~CUDAIntegrator() {

}

template<int dims>
void CUDAIntegrator<dims>::set_initial_rho(RhoMatrix<double> &r) {
    Integrator<dims>::set_initial_rho(r);
    _CPU_GPU();
}

template<int dims>
void CUDAIntegrator<dims>::_CPU_GPU() {
	for(unsigned int idx = 0; idx < this->_rho.bins(); idx++) {
		for(int species = 0; species < this->_model->N_species(); species++) {
			_h_rho(idx, species) = this->_rho(idx, species);
		}
	}

	CUDA_SAFE_CALL(cudaMemcpy(_d_rho, _h_rho.data(), _d_vec_size, cudaMemcpyHostToDevice));
}

template<int dims>
RhoMatrix<double> &CUDAIntegrator<dims>::rho() {
    if(!_output_ready) {
        CUDA_SAFE_CALL(cudaMemcpy(_h_rho.data(), _d_rho, _d_vec_size, cudaMemcpyDeviceToHost));

        for(unsigned int idx = 0; idx < this->_rho.bins(); idx++) {
            for(int species = 0; species < this->_model->N_species(); species++) {
                this->_rho(idx, species) = _h_rho(idx, species);
            }
        }

        _output_ready = true;
    }

    return this->_rho;
}

template class CUDAIntegrator<1>;
template class CUDAIntegrator<2>;

} /* namespace ch */
