#include "SemiImplicitMobilityCPU.h"

#include "../utils/utility_functions.h"

#include <cmath>
#include <random>

namespace ch {

template<int dims>
SemiImplicitMobilityCPU<dims>::SemiImplicitMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) :
        EulerCPU<dims>(sim_state, model, config) {

    _S = this->template _config_optional_value<double>(config, "semi_implicit.stabilization", 0.0);
    _rho_floor = this->template _config_optional_value<double>(config, "semi_implicit.rho_floor", 0.0);
    _cg_tol = this->template _config_optional_value<double>(config, "semi_implicit.cg_tol", 1e-8);
    _cg_max_iter = this->template _config_optional_value<int>(config, "semi_implicit.cg_max_iter", 500);

    _with_noise = this->template _config_optional_value<bool>(config, "mobility.with_noise", false);
    if(_with_noise) {
        double noise_rescale_factor = this->template _config_optional_value<double>(config, "mobility.noise_rescale_factor", 1.0);
        _noise_factor = std::sqrt(2.0 / (pow_dims<dims>(this->_dx) * this->_dt)) * noise_rescale_factor;
        this->info("Semi-implicit CH with non-constant mobility and noise (noise_factor = {})", _noise_factor);

        long long int seed = this->template _config_optional_value<long long int>(config, "seed", std::time(NULL));
        _generator.seed(seed);
    }
    else {
        this->info("Semi-implicit CH with non-constant mobility (no noise)");
    }

    if(_S > 0.0) {
        this->info("Semi-implicit stabilization enabled: S = {}", _S);
    }
    this->info("CG settings: tol = {}, max_iter = {}", _cg_tol, _cg_max_iter);
}

template<int dims>
SemiImplicitMobilityCPU<dims>::~SemiImplicitMobilityCPU() {
}

template<int dims>
double SemiImplicitMobilityCPU<dims>::_dot(const std::vector<double> &a, const std::vector<double> &b) {
    double s = 0.0;
    const size_t n = a.size();
    for(size_t i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

// Discrete laplacian for scalar arrays (periodic BC). We match EulerCPU's stencil.
template<>
void SemiImplicitMobilityCPU<1>::_laplacian_scalar(const std::vector<double> &in, std::vector<double> &out) {
    const int mask = this->_N_per_dim_minus_one;
    const unsigned int N = this->_N_bins;
    out.resize(N);

    for(unsigned int idx = 0; idx < N; idx++) {
        int idx_m = (int(idx) - 1 + int(N)) & mask;
        int idx_p = (int(idx) + 1) & mask;
        out[idx] = (in[idx_m] + in[idx_p] - 2.0 * in[idx]) / SQR(this->_dx);
    }
}

template<>
void SemiImplicitMobilityCPU<2>::_laplacian_scalar(const std::vector<double> &in, std::vector<double> &out) {
    const int mask = this->_N_per_dim_minus_one;
    const unsigned int N = this->_N_bins;
    out.resize(N);

    int coords[2];
    for(unsigned int idx = 0; idx < N; idx++) {
        this->_fill_coords(coords, (int)idx);

        int coords_xmy[2] = { (coords[0] - 1 + int(N)) & mask, coords[1] };
        int coords_xpy[2] = { (coords[0] + 1) & mask, coords[1] };
        int coords_xym[2] = { coords[0], (coords[1] - 1 + int(N)) & mask };
        int coords_xyp[2] = { coords[0], (coords[1] + 1) & mask };

        out[idx] = (
            in[this->_cell_idx(coords_xmy)] +
            in[this->_cell_idx(coords_xpy)] +
            in[this->_cell_idx(coords_xym)] +
            in[this->_cell_idx(coords_xyp)] -
            4.0 * in[idx]
        ) / SQR(this->_dx);
    }
}

// Build div( M grad(phi) ) into out, using a properly dimension-aware staggered discretization.
template<int dims>
void SemiImplicitMobilityCPU<dims>::_div_M_grad_from_scalar(const std::vector<double> &phi, int species, std::vector<double> &out) {
    static MultiField<Gradient<dims>> flux(this->_rho.bins(), this->_model->N_species());

    const int mask = this->_N_per_dim_minus_one;
    const unsigned int N = this->_N_bins;
    out.resize(N);

    int coords[dims];
    int coords_p[dims];

    for(unsigned int idx = 0; idx < N; idx++) {
        this->_fill_coords(coords, (int)idx);
        for(int d = 0; d < dims; d++) {
            memcpy(coords_p, coords, sizeof(int) * dims);
            coords_p[d] = (coords[d] + 1) & mask;
            int idx_p = this->_cell_idx(coords_p);

            const double M_idx = this->_sim_state.mobility(idx, species);
            const double M_p = this->_sim_state.mobility(idx_p, species);
            const double M_flux = 0.5 * (M_idx + M_p);

            flux(idx, species)[d] = M_flux * (phi[idx_p] - phi[idx]) / this->_dx;
        }
    }

    for(unsigned int idx = 0; idx < N; idx++) {
        out[idx] = this->_divergence(flux, species, idx);
    }
}

template<int dims>
void SemiImplicitMobilityCPU<dims>::_apply_D1(const std::vector<double> &x, int species, std::vector<double> &out) {
    _div_M_grad_from_scalar(x, species, out);
}

template<int dims>
void SemiImplicitMobilityCPU<dims>::_apply_D2(const std::vector<double> &x, int species, std::vector<double> &out) {
    static std::vector<double> lap;
    _laplacian_scalar(x, lap);
    _div_M_grad_from_scalar(lap, species, out);
}

template<int dims>
void SemiImplicitMobilityCPU<dims>::_apply_A(const std::vector<double> &x, int species, std::vector<double> &out) {
    static std::vector<double> D1;
    static std::vector<double> D2;

    _apply_D1(x, species, D1);
    _apply_D2(x, species, D2);

    const unsigned int N = this->_N_bins;
    out.resize(N);

    const double dt = this->_dt;
    const double k = this->_k_laplacian;

    for(unsigned int i = 0; i < N; i++) {
        out[i] = x[i] - dt * (_S * D1[i] - 2.0 * k * D2[i]);
    }
}

template<int dims>
int SemiImplicitMobilityCPU<dims>::_solve_CG(int species, const std::vector<double> &b, std::vector<double> &x) {
    const unsigned int N = this->_N_bins;
    x.resize(N);

    static std::vector<double> r;
    static std::vector<double> p;
    static std::vector<double> Ap;

    r.resize(N);
    p.resize(N);
    Ap.resize(N);

    // r = b - A x
    _apply_A(x, species, Ap);
    for(unsigned int i = 0; i < N; i++) {
        r[i] = b[i] - Ap[i];
        p[i] = r[i];
    }

    double rsold = _dot(r, r);
    const double rs0 = rsold;
    const double tol2 = _cg_tol * _cg_tol;

    if(rsold <= tol2 * (rs0 > 0.0 ? rs0 : 1.0)) return 0;

    for(int it = 0; it < _cg_max_iter; it++) {
        _apply_A(p, species, Ap);
        const double denom = _dot(p, Ap);
        if(std::fabs(denom) < 1e-30) return -1;

        const double alpha = rsold / denom;
        for(unsigned int i = 0; i < N; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        const double rsnew = _dot(r, r);
        if(rsnew <= tol2 * (rs0 > 0.0 ? rs0 : 1.0)) return it + 1;

        const double beta = rsnew / rsold;
        for(unsigned int i = 0; i < N; i++) p[i] = r[i] + beta * p[i];

        rsold = rsnew;
    }

    return -1;
}

template<int dims>
void SemiImplicitMobilityCPU<dims>::evolve() {
    static MultiField<double> rho_der(this->_rho.bins(), this->_model->N_species());
    static MultiField<double> rho_safe(this->_rho.bins(), this->_model->N_species());
    static MultiField<Gradient<dims>> stochastic_flux(this->_rho.bins(), this->_model->N_species());
    static std::normal_distribution<double> normal_dist(0.0, 1.0);

    // Optionally clamp densities before computing df_bulk/drho to avoid log instabilities
    if(_rho_floor > 0.0) {
        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            for(int s = 0; s < this->_model->N_species(); s++) {
                rho_safe(idx, s) = std::max(this->_rho(idx, s), _rho_floor);
            }
        }
        this->_model->der_bulk_free_energy(rho_safe, rho_der);
    }
    else {
        this->_model->der_bulk_free_energy(this->_rho, rho_der);
    }

    // Explicit chemical potential part: mu_exp = df_bulk/drho(rho^n) - S*rho^n
    static std::vector<double> mu_exp;
    static std::vector<double> rhs_div;
    static std::vector<double> noise_div;
    static std::vector<double> b;
    static std::vector<double> x;

    mu_exp.resize(this->_N_bins);
    rhs_div.resize(this->_N_bins);
    noise_div.assign(this->_N_bins, 0.0);
    b.resize(this->_N_bins);
    x.resize(this->_N_bins);

    // Optional noise contribution (kept explicit, as in EulerMobilityCPU)
    if(_with_noise) {
        const int mask = this->_N_per_dim_minus_one;
        int coords[dims], coords_p[dims];

        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            this->_fill_coords(coords, (int)idx);
            for(int species = 0; species < this->_model->N_species(); species++) {
                for(int d = 0; d < dims; d++) {
                    memcpy(coords_p, coords, sizeof(int) * dims);
                    coords_p[d] = (coords[d] + 1) & mask;
                    int idx_p = this->_cell_idx(coords_p);

                    const double M_idx = this->_sim_state.mobility(idx, species);
                    const double M_p = this->_sim_state.mobility(idx_p, species);
                    const double M_flux = 0.5 * (M_idx + M_p);

                    const double noise_amplitude = std::sqrt(M_flux) * _noise_factor;
                    stochastic_flux(idx, species)[d] = noise_amplitude * normal_dist(_generator);
                }
            }
        }

        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            for(int species = 0; species < this->_model->N_species(); species++) {
                // store per-species later (we overwrite noise_div per species below)
            }
        }
    }

    // Solve per species independently
    for(int species = 0; species < this->_model->N_species(); species++) {
        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            mu_exp[idx] = rho_der(idx, species) - _S * this->_rho(idx, species);
        }

        _div_M_grad_from_scalar(mu_exp, species, rhs_div);

        if(_with_noise) {
            for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
                noise_div[idx] = this->_divergence(stochastic_flux, species, idx);
            }
        }
        else {
            std::fill(noise_div.begin(), noise_div.end(), 0.0);
        }

        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            b[idx] = this->_rho(idx, species) + this->_dt * (rhs_div[idx] + noise_div[idx]);
            x[idx] = this->_rho(idx, species); // initial guess
        }

        const int iters = _solve_CG(species, b, x);
        if(iters < 0) this->warning("CG did not converge for species {}", species);

        for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
            this->_rho(idx, species) = x[idx];
        }
    }
}

template class SemiImplicitMobilityCPU<1>;
template class SemiImplicitMobilityCPU<2>;

} /* namespace ch */
