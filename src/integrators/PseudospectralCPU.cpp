/*
 * PseudospectralCPU.cpp
 *
 * Created on: 3/27/2024
 *     Author: Lorenzo
*/

#include "PseudospectralCPU.h"

namespace ch {

template<int dims>
PseudospectralCPU<dims>::PseudospectralCPU(FreeEnergyModel *model, toml::table &config) : Integrator<dims>(model, config) {
    _reciprocal_n.fill(this->_N_per_dim);
    hat_grid_size = _reciprocal_n[dims - 1] / 2 + 1; 
    for(int i = 0; i < dims - 1; i++) {
        hat_grid_size *= _reciprocal_n[i];
    }
    hat_vector_size = hat_grid_size * model->N_species();

    this->info("Size of the reciprocal vectors: {}", hat_vector_size);

    rho_hat.resize(hat_vector_size);
    rho_hat_copy.resize(hat_vector_size);
    sqr_wave_vectors.resize(hat_vector_size);
    dealiaser.resize(hat_vector_size);

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
            for(int kx_idx = 0; kx_idx < _reciprocal_n[0]; kx_idx++) {
                int kx = (kx_idx < _reciprocal_n[0] / 2) ? kx_idx : _reciprocal_n[0] - kx_idx;
                for(int ky = 0; ky < (_reciprocal_n[1] / 2 + 1); ky++) {
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
}

template<int dims>
PseudospectralCPU<dims>::~PseudospectralCPU() {
    fftw_destroy_plan(rho_inverse_plan);
    fftw_destroy_plan(f_der_plan);
    fftw_cleanup();
}

template<int dims>
void PseudospectralCPU<dims>::set_initial_rho(RhoMatrix<double> &r) {
    Integrator<dims>::set_initial_rho(r);

    f_der = RhoMatrix<double>(this->_rho.bins(), this->_model->N_species());
    f_der_hat.resize(hat_vector_size);
    
    int idist = this->_N_bins;
    int odist = hat_grid_size; // the distance between the first elements of arrays referring to nearby species in the reciprocal space
    f_der_plan = fftw_plan_many_dft_r2c(dims, _reciprocal_n.data(), this->_model->N_species(), f_der.data(), NULL, 1, idist, reinterpret_cast<fftw_complex *>(f_der_hat.data()), NULL, 1, odist, FFTW_ESTIMATE);

    // c2r transforms overwrite the input array
    rho_inverse_plan = fftw_plan_many_dft_c2r(dims, _reciprocal_n.data(), this->_model->N_species(), reinterpret_cast<fftw_complex *>(rho_hat_copy.data()), NULL, 1, odist, this->_rho.data(), NULL, 1, idist, FFTW_ESTIMATE);

    fftw_plan rho_plan = fftw_plan_many_dft_r2c(dims, _reciprocal_n.data(), this->_model->N_species(), this->_rho.data(), NULL, 1, idist, reinterpret_cast<fftw_complex *>(rho_hat.data()), NULL, 1, odist, FFTW_ESTIMATE);
    fftw_execute(rho_plan);
    rho_hat_copy = rho_hat;
    fftw_destroy_plan(rho_plan);
}

template<int dims>
void PseudospectralCPU<dims>::evolve() {
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_model->N_species(); species++) {
            f_der(idx, species) = this->_model->der_bulk_free_energy(species, this->_rho.rho_species(idx));
        }
    }

    fftw_execute(f_der_plan); // transform f_der into f_der_hat

    for(unsigned int k_idx = 0; k_idx < rho_hat.size(); k_idx++) {
        // f_der_hat[k_idx] *= dealiaser[k_idx];
        rho_hat[k_idx] = rho_hat_copy[k_idx] = (rho_hat[k_idx] - this->_dt * this->_M * sqr_wave_vectors[k_idx] * f_der_hat[k_idx]) / (1.0 + this->_dt * this->_M * 2.0 * this->_k_laplacian * SQR(sqr_wave_vectors[k_idx]));
        rho_hat_copy[k_idx] /= this->_N_bins;
    }

    fftw_execute(rho_inverse_plan);
}

template class PseudospectralCPU<1>;
template class PseudospectralCPU<2>;

} /* namespace ch */
