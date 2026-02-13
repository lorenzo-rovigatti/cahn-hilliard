/*
 * PseudospectralMobilityCPU.cpp
 *
 * Created on: 12/30/2025
 *     Author: Lorenzo
*/

#include "PseudospectralMobilityCPU.h"

namespace ch {

template<int dims>
PseudospectralMobilityCPU<dims>::PseudospectralMobilityCPU(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config) : 
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
    else if constexpr (dims == 3) {
        int k_idx = 0;
        for(int species = 0; species < model->N_species(); species++) {
            for(int kx_idx = 0; kx_idx < _reciprocal_n[0]; kx_idx++) {
                int kx_i = (kx_idx < _reciprocal_n[0] / 2) ? kx_idx : (kx_idx - _reciprocal_n[0]);
                double kx = 2.0 * M_PI * kx_i / (this->_N_per_dim * this->_dx);

                for(int ky_idx = 0; ky_idx < _reciprocal_n[1]; ky_idx++) {
                    int ky_i = (ky_idx < _reciprocal_n[1] / 2) ? ky_idx : (ky_idx - _reciprocal_n[1]);
                    double ky = 2.0 * M_PI * ky_i / (this->_N_per_dim * this->_dx);

                    for(int kz_idx = 0; kz_idx < (_reciprocal_n[2] / 2 + 1); kz_idx++) {
                        double kz = 2.0 * M_PI * kz_idx / (this->_N_per_dim * this->_dx);

                        double k2 = SQR(kx) + SQR(ky) + SQR(kz);
                        double k = std::sqrt(k2);
                        dealiaser[k_idx] = (k <= k_cut) ? 1.0 : 0.0;

                        kcomp[0][k_idx] = kx;
                        kcomp[1][k_idx] = ky;
                        kcomp[2][k_idx] = kz;
                        sqr_wave_vectors[k_idx] = k2;
                        k_idx++;
                    }
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

    // --- Option C workspaces ---
    tmp_real = MultiField<double>(this->_N_bins, this->_N_species);
    tmp_hat.resize(_hat_vector_size);
    q_hat.resize(_hat_vector_size);
    div_hat.resize(_hat_vector_size);

    M0_species.resize(this->_N_species, 0.0);

    // tmp plans: generic many r2c/c2r for tmp_real/tmp_hat
    tmp_r2c_plan = fftw_plan_many_dft_r2c(
        dims, _reciprocal_n.data(),
        this->_N_species,
        tmp_real.data(), nullptr, 1, idist,
        reinterpret_cast<fftw_complex*>(tmp_hat.data()), nullptr, 1, odist,
        FFTW_ESTIMATE
    );

    tmp_c2r_plan = fftw_plan_many_dft_c2r(
        dims, _reciprocal_n.data(),
        this->_N_species,
        reinterpret_cast<fftw_complex*>(tmp_hat.data()), nullptr, 1, odist,
        tmp_real.data(), nullptr, 1, idist,
        FFTW_ESTIMATE
    );

    // GMRES params from config (optional)
    gmres_restart = this->template _config_optional_value<int>(config, "pseudospectral.gmres_restart", 30);
    gmres_max_iter = this->template _config_optional_value<int>(config, "pseudospectral.gmres_max_iter", 200);
    gmres_tol = this->template _config_optional_value<double>(config, "pseudospectral.gmres_tol", 1e-10);

    _use_gmres = this->template _config_optional_value<bool>(config, "pseudospectral.use_gmres", false);
}

template<int dims>
PseudospectralMobilityCPU<dims>::~PseudospectralMobilityCPU() {
    fftw_destroy_plan(rho_inverse_plan);
    fftw_destroy_plan(f_der_plan);
    for(int d = 0; d < dims; d++) {
        fftw_destroy_plan(grad_mu_inverse_plan[d]);
        fftw_destroy_plan(flux_plan[d]);
    }

    fftw_destroy_plan(tmp_r2c_plan);
    fftw_destroy_plan(tmp_c2r_plan);

    fftw_cleanup();
}

template<int dims>
void PseudospectralMobilityCPU<dims>::evolve() {
    if(!_use_gmres) {
        _evolve_simple();
    } 
    else {
        _evolve_full();
    }
}

template<int dims>
void PseudospectralMobilityCPU<dims>::_evolve_simple() {
    // calculate reference mobility M0 per species
    _compute_M0_species();

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
            double dM = this->_sim_state.mobility(i, s) - M0_species[s];
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
        int s = k / _hat_grid_size;
        double k2 = sqr_wave_vectors[k];
        double k4 = SQR(k2);

        double denom = 1.0 + this->_dt * M0_species[s] * (_S * k2 + 2.0 * this->_k_laplacian * k4);
        rho_hat[k] = rho_hat_copy[k] = (rho_hat[k] + this->_dt * divJ_hat[k] - this->_dt * M0_species[s] * k2 * (f_der_hat[k] - _S * rho_hat[k])) / denom;

        rho_hat_copy[k] /= this->_N_bins;
    }

    // 8) back to real space
    fftw_execute(rho_inverse_plan);
}

template<int dims>
void PseudospectralMobilityCPU<dims>::_evolve_full() {
    // 0) compute per-species M0 for preconditioner (max over space)
    _compute_M0_species();

    // 1) mu_bulk(rho^n) in real space into f_der
    this->_model->der_bulk_free_energy(this->_rho, f_der);

    // 2) Build g = mu_bulk - S*rho^n in tmp_real (real space)
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            tmp_real(i, s) = f_der(i, s) - _S * this->_rho(i, s);
        }
    }

    // 3) Compute rhs b = rho^n + dt * div( M grad g )
    //    We'll compute div(M grad g) into tmp_real (overwrite), then build b_vec.
    //    Steps: FFT g -> tmp_hat, grad in real via grad_mu buffers, flux = M*grad, div in tmp_real.
    fftw_execute(tmp_r2c_plan); // tmp_real (g) -> tmp_hat

    // grad g in real: grad_mu[d] (reuse buffers)
    for(int d = 0; d < dims; d++) {
        for(int k = 0; k < _hat_vector_size; k++) {
            grad_mu_hat[d][k] = std::complex<double>(0.0, kcomp[d][k]) * tmp_hat[k];
            if(use_dealias) grad_mu_hat[d][k] *= dealiaser[k];
        }
        fftw_execute(grad_mu_inverse_plan[d]);

        for(int i = 0; i < this->_N_bins; i++)
            for(int s = 0; s < this->_N_species; s++)
                grad_mu[d](i, s) /= this->_N_bins; // normalize c2r
    }

    // flux = M * grad g
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            double Mxs = this->_sim_state.mobility(i, s);
            for(int d = 0; d < dims; d++)
                flux[d](i, s) = Mxs * grad_mu[d](i, s);
        }
    }

    // div_hat = i k · flux_hat
    std::fill(div_hat.begin(), div_hat.end(), 0.0);
    for(int d = 0; d < dims; d++) {
        fftw_execute(flux_plan[d]); // flux[d] -> flux_hat[d]
        for(int k = 0; k < _hat_vector_size; k++) {
            if(use_dealias) flux_hat[d][k] *= dealiaser[k];
            div_hat[k] += std::complex<double>(0.0, kcomp[d][k]) * flux_hat[d][k];
        }
    }

    // inverse FFT div_hat -> tmp_real (div in real)
    tmp_hat = div_hat;
    fftw_execute(tmp_c2r_plan);
    for(int i = 0; i < this->_N_bins; i++)
        for(int s = 0; s < this->_N_species; s++)
            tmp_real(i, s) /= this->_N_bins;

    // build b vector
    const int n = this->_N_bins * this->_N_species;
    std::vector<double> b(n), x(n);

    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            int idx = s * this->_N_bins + i;
            b[idx] = this->_rho(i, s) + this->_dt * tmp_real(i, s);
            x[idx] = this->_rho(i, s); // initial guess: previous rho
        }
    }

    // 4) Solve A x = b with left-preconditioned restarted GMRES
    int iters = _gmres_solve(b, x);
    (void)iters; // you can log it if you want

    // 5) Copy solution back to rho (real space)
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            this->_rho(i, s) = x[s * this->_N_bins + i];
        }
    }

    // 6) Update rho_hat from rho for next step
    //    Use tmp_real/tmp_hat plans: copy rho -> tmp_real, FFT -> tmp_hat, then rho_hat = tmp_hat
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            tmp_real(i, s) = this->_rho(i, s);
        }
    }

    fftw_execute(tmp_r2c_plan);
    rho_hat = tmp_hat;
    rho_hat_copy = rho_hat;
}

template<int dims>
void PseudospectralMobilityCPU<dims>::_compute_M0_species() {
    std::fill(M0_species.begin(), M0_species.end(), 0.0);
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            M0_species[s] = std::max(M0_species[s], this->_sim_state.mobility(i, s));
        }
    }
}

template<int dims>
double PseudospectralMobilityCPU<dims>::_dot(const std::vector<double>& a, const std::vector<double>& b) const {
    double sum = 0.0;
    for(size_t i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

template<int dims>
double PseudospectralMobilityCPU<dims>::_norm2(const std::vector<double>& a) const {
    return std::sqrt(_dot(a, a));
}

template<int dims>
void PseudospectralMobilityCPU<dims>::_apply_Pinv(const std::vector<double> &r, std::vector<double> &z) {
    z.resize(r.size());

    // copy r -> tmp_real
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            tmp_real(i, s) = r[s * this->_N_bins + i];
        }
    }

    // FFT
    fftw_execute(tmp_r2c_plan); // tmp_real -> tmp_hat

    // divide by diagonal in k-space
    for(int k = 0; k < _hat_vector_size; k++) {
        int s = k / _hat_grid_size;
        double k2v = sqr_wave_vectors[k];
        double k4v = k2v * k2v;

        double denom = 1.0 + this->_dt * M0_species[s] * (_S * k2v + 2.0 * this->_k_laplacian * k4v);
        tmp_hat[k] /= denom;
    }

    // inverse FFT
    fftw_execute(tmp_c2r_plan); // tmp_hat -> tmp_real (unnormalized)
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            z[s * this->_N_bins + i] = tmp_real(i, s) / this->_N_bins;
        }
    }
}

template<int dims>
void PseudospectralMobilityCPU<dims>::_apply_A(const std::vector<double> &x, std::vector<double> &Ax) {
    Ax.resize(x.size());

    // Copy x -> tmp_real
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            tmp_real(i, s) = x[s * this->_N_bins + i];
        }
    }

    // FFT x -> tmp_hat
    fftw_execute(tmp_r2c_plan);

    // q_hat = (S + 2k*k^2) * x_hat   (since -2k laplacian x => +2k*k^2 * x_hat)
    for(int k = 0; k < _hat_vector_size; k++) {
        double k2v = sqr_wave_vectors[k];
        q_hat[k] = (_S + 2.0 * this->_k_laplacian * k2v) * tmp_hat[k];
        if(use_dealias) {
            q_hat[k] *= dealiaser[k];
        }
    }

    // grad q in real
    for(int d = 0; d < dims; d++) {
        for(int k = 0; k < _hat_vector_size; k++) {
            grad_mu_hat[d][k] = std::complex<double>(0.0, kcomp[d][k]) * q_hat[k];
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

    // flux = M * grad q
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            double Mxs = this->_sim_state.mobility(i, s);
            for(int d = 0; d < dims; d++) {
                flux[d](i, s) = Mxs * grad_mu[d](i, s);
            }
        }
    }

    // div_hat = i k · flux_hat
    std::fill(div_hat.begin(), div_hat.end(), 0.0);
    for(int d = 0; d < dims; d++) {
        fftw_execute(flux_plan[d]);
        for(int k = 0; k < _hat_vector_size; k++) {
            if(use_dealias) flux_hat[d][k] *= dealiaser[k];
            div_hat[k] += std::complex<double>(0.0, kcomp[d][k]) * flux_hat[d][k];
        }
    }

    // div_hat -> tmp_real
    tmp_hat = div_hat;
    fftw_execute(tmp_c2r_plan);

    // Ax = x - dt * div (normalize c2r)
    for(int i = 0; i < this->_N_bins; i++) {
        for(int s = 0; s < this->_N_species; s++) {
            double div_real = tmp_real(i, s) / this->_N_bins;
            Ax[s * this->_N_bins + i] = x[s * this->_N_bins + i] - this->_dt * div_real;
        }
    }
}

template<int dims>
int PseudospectralMobilityCPU<dims>::_gmres_solve(const std::vector<double> &b,
                                                 std::vector<double> &x) {
    const int n = (int)b.size();
    const int m = std::max(2, gmres_restart);

    // rhs = P^{-1} b
    std::vector<double> rhs;
    _apply_Pinv(b, rhs);

    auto apply_M = [&](const std::vector<double>& v, std::vector<double>& Mv) {
        std::vector<double> Av, tmp;
        _apply_A(v, Av);
        _apply_Pinv(Av, Mv); // left preconditioning
    };

    // r0 = rhs - M x0
    std::vector<double> Mx, r;
    apply_M(x, Mx);
    r.resize(n);
    for(int i = 0; i < n; i++) r[i] = rhs[i] - Mx[i];

    const double rhs_norm = std::max(_norm2(rhs), 1e-30);
    double beta = _norm2(r);
    double rel = beta / rhs_norm;
    if(rel < gmres_tol) {
        return 0;
    }

    int iter_total = 0;

    // Restart loop
    while(iter_total < gmres_max_iter) {
        // V basis vectors (m+1)
        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n, 0.0));
        // Hessenberg
        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
        // Givens
        std::vector<double> cs(m, 0.0), sn(m, 0.0);
        // residual in Krylov coords
        std::vector<double> g(m + 1, 0.0);

        // V0 = r / beta
        for(int i = 0; i < n; i++) V[0][i] = r[i] / beta;
        g[0] = beta;

        int j_end = 0;

        for(int j = 0; j < m && iter_total < gmres_max_iter; j++, iter_total++) {
            // w = M V[j]
            std::vector<double> w;
            apply_M(V[j], w);

            // Arnoldi
            for(int i = 0; i <= j; i++) {
                H[i][j] = _dot(w, V[i]);
                for(int k = 0; k < n; k++) w[k] -= H[i][j] * V[i][k];
            }
            H[j + 1][j] = _norm2(w);

            if(H[j + 1][j] > 1e-30) {
                for(int k = 0; k < n; k++) V[j + 1][k] = w[k] / H[j + 1][j];
            } else {
                // happy breakdown
                j_end = j + 1;
                break;
            }

            // Apply previous Givens rotations
            for(int i = 0; i < j; i++) {
                double tmp1 = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                double tmp2 = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = tmp1;
                H[i + 1][j] = tmp2;
            }

            // Create new Givens to zero H[j+1][j]
            double a = H[j][j];
            double b2 = H[j + 1][j];
            double denom = std::sqrt(a*a + b2*b2);
            cs[j] = (denom > 0) ? (a / denom) : 1.0;
            sn[j] = (denom > 0) ? (b2 / denom) : 0.0;

            // Apply it
            H[j][j] = cs[j] * a + sn[j] * b2;
            H[j + 1][j] = 0.0;

            // Update g
            double gj = g[j];
            g[j]     = cs[j] * gj;
            g[j + 1] = -sn[j] * gj;

            rel = std::abs(g[j + 1]) / rhs_norm;
            if(rel < gmres_tol) {
                j_end = j + 1;
                break;
            }

            j_end = j + 1;
        }

        // Solve y from upper triangular system Hy = g (size j_end)
        std::vector<double> y(j_end, 0.0);
        for(int i = j_end - 1; i >= 0; i--) {
            double sum = g[i];
            for(int k = i + 1; k < j_end; k++) sum -= H[i][k] * y[k];
            y[i] = sum / H[i][i];
        }

        // Update x = x + V(:,0..j_end-1) * y
        for(int i = 0; i < j_end; i++) {
            for(int k = 0; k < n; k++) x[k] += V[i][k] * y[i];
        }

        // Recompute residual r = rhs - Mx
        apply_M(x, Mx);
        for(int i = 0; i < n; i++) r[i] = rhs[i] - Mx[i];
        beta = _norm2(r);
        rel = beta / rhs_norm;

        if(rel < gmres_tol) break;
    }

    return iter_total;
}

template class PseudospectralMobilityCPU<1>;
template class PseudospectralMobilityCPU<2>;
template class PseudospectralMobilityCPU<3>;

} /* namespace ch */