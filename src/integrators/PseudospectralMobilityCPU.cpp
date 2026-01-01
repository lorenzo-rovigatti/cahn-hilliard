/*
 * PseudospectralMobilityCPU.cpp
 *
 * Created on: 12/30/2025
 *     Author: Lorenzo
*/

#include "PseudospectralMobilityCPU.h"

namespace ch {

template<int dims>
PseudospectralMobilityCPU<dims>::PseudospectralMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
        Integrator<dims>(sim_state, model, config) {
    _reciprocal_n.fill(this->_N_per_dim);
    _hat_grid_size = _reciprocal_n[dims - 1] / 2 + 1;
    for(int i = 0; i < dims - 1; i++) {
        _hat_grid_size *= _reciprocal_n[i];
    }
    _hat_vector_size = _hat_grid_size * model->N_species();

    _S = this->template _config_optional_value<double>(config, "pseudospectral.S", 0.0);
    use_dealias = this->template _config_optional_value<bool>(config, "pseudospectral.use_dealias", false);

    this->info("Size of the reciprocal vectors: {}, S = {}, use_dealias = {}", _hat_vector_size, _S, use_dealias);

    rho_hat.resize(_hat_vector_size);
    rho_hat_copy.resize(_hat_vector_size);

    // bulk mu = df/d(rho)
    f_der = MultiField<double>(this->_rho.bins(), this->_model->N_species());
    f_der_hat.resize(_hat_vector_size);

    // mu_hat buffer
    mu_hat.resize(_hat_vector_size);

    // divergence buffer (hat)
    divJ_hat.resize(_hat_vector_size);

    // dealias / k2 / k-components
    dealiaser.resize(_hat_vector_size);
    sqr_wave_vectors.resize(_hat_vector_size);
    for(int d = 0; d < dims; d++) {
        kcomp[d].resize(_hat_vector_size);
    }

    // real-space grad mu and flux (one field per dimension)
    for(int d = 0; d < dims; d++) {
        grad_mu[d] = MultiField<double>(this->_rho.bins(), this->_model->N_species());
        flux[d]    = MultiField<double>(this->_rho.bins(), this->_model->N_species());
        flux_hat[d].resize(_hat_vector_size);
        grad_mu_hat[d].resize(_hat_vector_size);
    }

    // Build wavevectors and dealiaser
    double k_cut = this->_N_per_dim * M_PI / (this->_N_per_dim * this->_dx) * 2.0 / 3.0;

    if constexpr (dims == 1) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(unsigned int i = 0; i < _hat_grid_size; i++) {
                double kx = 2.0 * M_PI * i / (this->_N_per_dim * this->_dx);
                double k = std::abs(kx);
                dealiaser[k_idx] = (k <= k_cut) ? 1.0 : 0.0;

                kcomp[0][k_idx] = kx;
                sqr_wave_vectors[k_idx] = SQR(kx);
                k_idx++;
            }
        }
    } 
    else if constexpr (dims == 2) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(int kx_idx = 0; kx_idx < _reciprocal_n[0]; kx_idx++) {
                int kx_i = (kx_idx < _reciprocal_n[0] / 2) ? kx_idx : (kx_idx - _reciprocal_n[0]);
                double kx = 2.0 * M_PI * kx_i / (this->_N_per_dim * this->_dx);

                for(int ky_idx = 0; ky_idx < (_reciprocal_n[1] / 2 + 1); ky_idx++) {
                    double ky = 2.0 * M_PI * ky_idx / (this->_N_per_dim * this->_dx);

                    double k2 = SQR(kx) + SQR(ky);
                    double k = std::sqrt(k2);
                    dealiaser[k_idx] = (k <= k_cut) ? 1.0 : 0.0;

                    kcomp[0][k_idx] = kx;
                    kcomp[1][k_idx] = ky;
                    sqr_wave_vectors[k_idx] = k2;
                    k_idx++;
                }
            }
        }
    } 
    else {
        this->critical("Unsupported number of dimensions {}", dims);
    }

    // FFTW plans layout
    int idist = this->_N_bins;
    int odist = _hat_grid_size;

    // df/drho real -> hat
    f_der_plan = fftw_plan_many_dft_r2c(
        dims, _reciprocal_n.data(),
        this->_model->N_species(),
        f_der.data(), nullptr, 1, idist,
        reinterpret_cast<fftw_complex*>(f_der_hat.data()), nullptr, 1, odist, FFTW_ESTIMATE
    );

    // rho hat -> rho real (note: c2r overwrites input, so we use rho_hat_copy)
    rho_inverse_plan = fftw_plan_many_dft_c2r(
        dims, _reciprocal_n.data(),
        this->_model->N_species(),
        reinterpret_cast<fftw_complex*>(rho_hat_copy.data()), nullptr, 1, odist,
        this->_rho.data(), nullptr, 1, idist, FFTW_ESTIMATE
    );

    // grad_mu_hat[d] (complex) -> grad_mu[d] (real)
    for(int d = 0; d < dims; d++) {
        grad_mu_inverse_plan[d] = fftw_plan_many_dft_c2r(
            dims, _reciprocal_n.data(),
            this->_model->N_species(),
            reinterpret_cast<fftw_complex*>(grad_mu_hat[d].data()), nullptr, 1, odist,
            grad_mu[d].data(), nullptr, 1, idist,
            FFTW_ESTIMATE
        );

        // flux[d] real -> flux_hat[d] complex
        flux_plan[d] = fftw_plan_many_dft_r2c(
            dims, _reciprocal_n.data(),
            this->_model->N_species(),
            flux[d].data(), nullptr, 1, idist,
            reinterpret_cast<fftw_complex*>(flux_hat[d].data()), nullptr, 1, odist,
            FFTW_ESTIMATE
        );
    }

    // Initialize rho_hat from current rho
    fftw_plan rho_plan = fftw_plan_many_dft_r2c(
        dims, _reciprocal_n.data(),
        this->_model->N_species(),
        this->_rho.data(), nullptr, 1, idist,
        reinterpret_cast<fftw_complex*>(rho_hat.data()), nullptr, 1, odist,
        FFTW_ESTIMATE
    );
    fftw_execute(rho_plan);
    rho_hat_copy = rho_hat;
    fftw_destroy_plan(rho_plan);
}

template<int dims>
PseudospectralMobilityCPU<dims>::~PseudospectralMobilityCPU() {
    fftw_destroy_plan(rho_inverse_plan);
    fftw_destroy_plan(f_der_plan);
    for(int d = 0; d < dims; d++) {
        fftw_destroy_plan(grad_mu_inverse_plan[d]);
        fftw_destroy_plan(flux_plan[d]);
    }
    fftw_cleanup();
}

template<int dims>
void PseudospectralMobilityCPU<dims>::evolve() {
    // calculate reference mobility M0 (max over space/species is safest)
    _M0 = 0.0;
    for(int idx = 0; idx < this->_N_bins; idx++) {
        for(int s = 0; s < this->_model->N_species(); s++) {
            _M0 = std::max(_M0, this->_sim_state.mobility(idx, s));
        }
    }

    // 1) bulk chemical potential in real space
    this->_model->der_bulk_free_energy(this->_rho, f_der);

    // 2) FFT(mu_bulk)
    fftw_execute(f_der_plan);

    // 3) build full mu_hat = mu_bulk_hat + 2k * k^2 * rho_hat
    for(int k = 0; k < _hat_vector_size; k++) {
        mu_hat[k] = f_der_hat[k] + 2.0 * this->_k_laplacian * sqr_wave_vectors[k] * rho_hat[k];
        if(use_dealias) {
            mu_hat[k] *= dealiaser[k];
        }
    }

    // 4) grad(mu) in real space
    for(int d = 0; d < dims; d++) {
        for(int k = 0; k < _hat_vector_size; k++) {
            grad_mu_hat[d][k] = std::complex<double>(0.0, kcomp[d][k]) * mu_hat[k];
            if(use_dealias) {
                grad_mu_hat[d][k] *= dealiaser[k];
            }
        }

        fftw_execute(grad_mu_inverse_plan[d]);

        for(int i = 0; i < this->_N_bins; i++) {
            for(int s = 0; s < this->_N_species; s++) {
                grad_mu[d](i, s) /= this->_N_bins;
            }
        }
    }

    // 5) real-space flux: J = (M - M0) grad(mu)
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            double dM = this->_sim_state.mobility(i, s) - _M0;
            for(int d = 0; d < dims; d++) {
                flux[d](i, s) = dM * grad_mu[d](i, s);
            }
        }
    }

    // 6) FFT and divergence
    std::fill(divJ_hat.begin(), divJ_hat.end(), 0.0);

    for(int d = 0; d < dims; d++) {
        fftw_execute(flux_plan[d]);
        for(int k = 0; k < _hat_vector_size; k++) {
            if(use_dealias) {
                flux_hat[d][k] *= dealiaser[k];
            }
            divJ_hat[k] += std::complex<double>(0.0, kcomp[d][k]) * flux_hat[d][k];
        }
    }

    // 7) stabilized semi-implicit update
    for(int k = 0; k < _hat_vector_size; k++) {
        double k2 = sqr_wave_vectors[k];
        double k4 = SQR(k2);

        double denom = 1.0 + this->_dt * _M0 * (_S * k2 + 2.0 * this->_k_laplacian * k4);

        rho_hat[k] = rho_hat_copy[k] = (rho_hat[k] + this->_dt * divJ_hat[k] - this->_dt * _M0 * k2 * (f_der_hat[k] - _S * rho_hat[k])) / denom;

        rho_hat_copy[k] /= this->_N_bins;
    }

    // 8) back to real space
    fftw_execute(rho_inverse_plan);
}


template class PseudospectralMobilityCPU<1>;
template class PseudospectralMobilityCPU<2>;

} /* namespace ch */