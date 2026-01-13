/*
 * PseudospectralCPU.cpp
 *
 * Created on: 3/27/2024
 *     Author: Lorenzo
*/

#include "PseudospectralCPU.h"

namespace ch {

template<int dims>
PseudospectralCPU<dims>::PseudospectralCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
        Integrator<dims>(sim_state, model, config) {

    _reciprocal_n.fill(this->_N_per_dim);
    hat_grid_size = _reciprocal_n[dims - 1] / 2 + 1; 
    for(int i = 0; i < dims - 1; i++) {
        hat_grid_size *= _reciprocal_n[i];
    }
    hat_vector_size = hat_grid_size * model->N_species();

    _splitting_S = this->template _config_optional_value<double>(config, "pseudospectral.S", 0.0);
    use_dealias = this->template _config_optional_value<bool>(config, "pseudospectral.use_dealias", false);

    this->info("Size of the reciprocal vectors: {}, S = {}, use_dealias = {}", hat_vector_size, _splitting_S, use_dealias);

    rho_hat.resize(hat_vector_size);
    rho_hat_copy.resize(hat_vector_size);
    sqr_wave_vectors.resize(hat_vector_size);
    dealiaser.resize(hat_vector_size);

    double k_cut = this->_N_per_dim * M_PI / (this->_N_per_dim * this->_dx) * 2.0 / 3.0;
    if constexpr (dims == 1) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(unsigned int i = 0; i < hat_grid_size; i++) {
                double k = 2.0 * M_PI * i / (this->_N_per_dim * this->_dx);
                dealiaser[k_idx] = (k <= k_cut) ? 1.0 : 0.0;
                sqr_wave_vectors[k_idx] = SQR(k);
                k_idx++;
            }
        }
    }
    else if constexpr (dims == 2) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(int kx_idx = 0; kx_idx < _reciprocal_n[0]; kx_idx++) {
                int kx = (kx_idx < _reciprocal_n[0] / 2) ? kx_idx : _reciprocal_n[0] - kx_idx;
                for(int ky = 0; ky < (_reciprocal_n[1] / 2 + 1); ky++) {
                    double k = 2.0 * M_PI * std::sqrt(SQR(kx) + SQR(ky)) / (this->_N_per_dim * this->_dx);
                    dealiaser[k_idx] = (k < k_cut) ? 1.0 : 0.0;
                    sqr_wave_vectors[k_idx] = SQR(k);
                    k_idx++;
                }
            }
        }
    }
    else {
        this->critical("Unsupported number of dimensions {}", dims);
    }

    f_der = MultiField<double>(this->_rho.bins(), this->_model->N_species());
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
PseudospectralCPU<dims>::~PseudospectralCPU() {
    fftw_destroy_plan(rho_inverse_plan);
    fftw_destroy_plan(f_der_plan);
    fftw_cleanup();
}

template<int dims>
void PseudospectralCPU<dims>::evolve() {
    static double M = this->_sim_state.mobility(0, 0); // constant mobility

    this->_model->der_bulk_free_energy(this->_rho, f_der);

    fftw_execute(f_der_plan); // transform f_der into f_der_hat

    for(unsigned int k_idx = 0; k_idx < rho_hat.size(); k_idx++) {
        if(use_dealias) {
            f_der_hat[k_idx] *= dealiaser[k_idx];
        }

        double k2 = sqr_wave_vectors[k_idx];
        double k4 = SQR(k2);

        double denom = 1.0 + this->_dt * M * (_splitting_S * k2 + 2.0 * this->_k_laplacian * k4);
        rho_hat[k_idx] = rho_hat_copy[k_idx] = (rho_hat[k_idx] - this->_dt * M * k2 * (f_der_hat[k_idx] - _splitting_S * rho_hat[k_idx])) / denom;

        rho_hat_copy[k_idx] /= this->_N_bins;

    }

    fftw_execute(rho_inverse_plan);
}

template class PseudospectralCPU<1>;
template class PseudospectralCPU<2>;

} /* namespace ch */
