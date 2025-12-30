#include "SemiImplicitMobilityCPU.h"

#include "../utils/Gradient.h"
#include "../utils/utility_functions.h"

#include <cmath>
#include <random>

namespace ch {

template<int dims>
SemiImplicitMobilityCPU<dims>::SemiImplicitMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) :
        PseudospectralCPU<dims>(sim_state, model, config) {
    _N_per_dim_minus_one = this->_N_per_dim - 1;
    // reference mobility
    _M0 = this->template _config_optional_value<double>(config, "mobility.M0", this->_sim_state.mobility(0,0));

    _rho_floor = this->template _config_optional_value<double>(config, "semi_implicit.rho_floor", 0.0);

    use_dealias = this->template _config_optional_value<bool>(config, "semi_implicit.dealias", false);

    mu_real = MultiField<double>(this->_rho.bins(), model->N_species());
    corr_real = MultiField<double>(this->_rho.bins(), model->N_species());
    corr_hat.resize(this->hat_vector_size);

    this->info("SemiImplicitMobilityCPU: M0={}, rho_floor={}, dealias={}", _M0, _rho_floor, use_dealias);
}

template<int dims>
SemiImplicitMobilityCPU<dims>::~SemiImplicitMobilityCPU() {}

/* -------- helpers: periodic indexing like EulerCPU -------- */
template<int dims>
void SemiImplicitMobilityCPU<dims>::_fill_coords(int coords[dims], int idx) const {
    // Works if _N_per_dim is power-of-two (as in your codebase).
    for(int d = 0; d < dims; d++) {
        coords[d] = idx & _N_per_dim_minus_one;
        idx >>= int(std::log2(this->_N_per_dim));
    }
}

template<int dims>
int SemiImplicitMobilityCPU<dims>::_cell_idx(const int coords[dims]) const {
    int idx = 0;
    for(int d = dims - 1; d >= 0; d--) {
        idx <<= int(std::log2(this->_N_per_dim));
        idx |= (coords[d] & _N_per_dim_minus_one);
    }
    return idx;
}

template<>
double SemiImplicitMobilityCPU<1>::_cell_laplacian(MultiField<double> &field, int species, int idx)  const {
	int idx_m = (idx - 1 + this->_N_bins) & _N_per_dim_minus_one;
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(this->_dx);
}

template<>
double SemiImplicitMobilityCPU<2>::_cell_laplacian(MultiField<double> &field, int species, int idx) const {
	int coords_xy[2];
	_fill_coords(coords_xy, idx);

	int coords_xmy[2] = {
			(coords_xy[0] - 1 + this->_N_bins) & _N_per_dim_minus_one,
			coords_xy[1]
	};

	int coords_xym[2] = {
			coords_xy[0],
			(coords_xy[1] - 1 + this->_N_bins) & _N_per_dim_minus_one
	};

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & _N_per_dim_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & _N_per_dim_minus_one
	};

	return (
			field(_cell_idx(coords_xmy), species) +
			field(_cell_idx(coords_xpy), species) +
			field(_cell_idx(coords_xym), species) +
			field(_cell_idx(coords_xyp), species) -
			4 * field(idx, species))
			/ SQR(this->_dx);
}

// ------------------------------------------------------------
// μ = f'(ρ) - 2κ ∇²ρ   (real space)
// ------------------------------------------------------------
template<int dims>
void SemiImplicitMobilityCPU<dims>::compute_mu() {

    // f'(rho) already stored in f_der by base class logic
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int s = 0; s < this->_model->N_species(); s++) {
            mu_real(idx, s) = this->f_der(idx, s)
                - 2.0 * this->_k_laplacian
                  * this->_cell_laplacian(this->_rho, s, idx);
        }
    }
}

// ------------------------------------------------------------
// corr = div( (M - M0) grad(mu) )   (real space, FD)
// ------------------------------------------------------------
template<int dims>
void SemiImplicitMobilityCPU<dims>::compute_correction(int species) {

    int coords[dims], coords_p[dims], coords_m[dims];

    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {

        this->_fill_coords(coords, idx);
        double div = 0.0;

        for(int d = 0; d < dims; d++) {

            std::copy(coords, coords + dims, coords_p);
            std::copy(coords, coords + dims, coords_m);

            coords_p[d] = (coords[d] + 1) & (this->_N_per_dim - 1);
            coords_m[d] = (coords[d] - 1 + this->_N_per_dim) & (this->_N_per_dim - 1);

            int idx_p = this->_cell_idx(coords_p);
            int idx_m = this->_cell_idx(coords_m);

            double M_p = this->_sim_state.mobility(idx_p, species);
            double M_m = this->_sim_state.mobility(idx_m, species);

            double grad_p = (mu_real(idx_p, species) - mu_real(idx, species)) / this->_dx;
            double grad_m = (mu_real(idx, species) - mu_real(idx_m, species)) / this->_dx;

            div += ( (M_p - _M0) * grad_p - (M_m - _M0) * grad_m ) / this->_dx;
        }

        corr_real(idx, species) = div;
    }
}


// ------------------------------------------------------------
// evolve()
// ------------------------------------------------------------
template<int dims>
void SemiImplicitMobilityCPU<dims>::evolve() {

    // ---------- f'(rho) ----------
    if(_rho_floor > 0.0) {
        static MultiField<double> rho_safe(this->_rho.bins(), this->_model->N_species());
        for(unsigned int i = 0; i < this->_N_bins; i++)
            for(int s = 0; s < this->_model->N_species(); s++)
                rho_safe(i,s) = std::max(this->_rho(i,s), _rho_floor);

        this->_model->der_bulk_free_energy(rho_safe, this->f_der);
    } else {
        this->_model->der_bulk_free_energy(this->_rho, this->f_der);
    }

    // FFT f_der → f_der_hat
    fftw_execute(this->f_der_plan);

    // ---------- μ in real space ----------
    compute_mu();

    // ---------- correction term ----------
    for(int s = 0; s < this->_model->N_species(); s++)
        compute_correction(s);

    // FFT corr_real → corr_hat (reuse f_der buffer & plan)
    for(unsigned int i = 0; i < this->_N_bins; i++)
        for(int s = 0; s < this->_model->N_species(); s++)
            this->f_der(i,s) = corr_real(i,s);

    fftw_execute(this->f_der_plan);
    corr_hat = this->f_der_hat;

    // ---------- spectral update ----------
    for(unsigned int k_idx = 0; k_idx < this->rho_hat.size(); k_idx++) {

        double k2 = this->sqr_wave_vectors[k_idx];
        double k4 = SQR(k2);

        std::complex<double> rhs =
              this->rho_hat[k_idx]
            - this->_dt * _M0 * k2 * this->f_der_hat[k_idx]
            + this->_dt * corr_hat[k_idx];

        double denom = 1.0 + this->_dt * _M0 * 2.0 * this->_k_laplacian * k4;

        this->rho_hat[k_idx] = this->rho_hat_copy[k_idx] = rhs / denom;

        if(use_dealias)
            this->rho_hat_copy[k_idx] *= this->dealiaser[k_idx];

        this->rho_hat_copy[k_idx] /= this->_N_bins;
    }

    // ---------- inverse FFT ----------
    fftw_execute(this->rho_inverse_plan);
}

template class SemiImplicitMobilityCPU<1>;
template class SemiImplicitMobilityCPU<2>;

} /* namespace ch */
