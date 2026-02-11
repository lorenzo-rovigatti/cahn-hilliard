/*
 * BailoFiniteVolume.cpp
 *
 * Created on: 4/1/2024
 *     Author: Lorenzo
*/

#include "BailoFiniteVolume.h"

#include <fsolve/fsolve.hpp>

namespace ch {

int N_species;
int N_bins;
double k_laplacian;
double dx;
double dt;
FreeEnergyModel *fe_model;
MultiField<double> rho_curr;
MultiField<double> base_csi;

template<int dims> BailoFiniteVolume<dims> *bailo;

double laplacian_1D(double *rho, int species, int idx) {
	int idx_m = (idx - 1 + N_bins) & (N_bins - 1);
	int idx_p = (idx + 1) & (N_bins - 1);

	return (rho[species * N_bins + idx_m] + rho[species * N_bins + idx_p] - 2.0 * rho[species * N_bins + idx]) / SQR(dx);
}

template<int dims>
void to_solve(int n, double *rho, double *fvec) {
	static MultiField<double> csi(N_bins, N_species);
    static MultiField<double> F_half(N_bins, N_species);

    for(int idx = 0; idx < N_bins; idx++) {
		SpeciesView<double> rho_species(rho + idx, csi.bins(), N_species);
        for(int species = 0; species < N_species; species++) {
			double F_der_con = fe_model->der_bulk_free_energy_contractive(species, rho_species);
			double interf_con = 2 * k_laplacian * laplacian_1D(rho, species, idx);
			csi(idx, species) = base_csi(idx, species) + F_der_con - 0.5 * interf_con;
        }
    }

	for(int idx = 0; idx < N_bins; idx++) {
		int prev_idx = (idx == 0) ? N_bins - 1 : idx - 1;
        for(int species = 0; species < N_species; species++) {
			F_half(idx, species) = -(csi(prev_idx, species) - csi(idx, species)) / dx;
		}
	}

	for(int idx = 0; idx < N_bins; idx++) {
		int next_idx = (idx == N_bins - 1) ?  0 : idx + 1;
        for(int species = 0; species < N_species; species++) {
			int matrix_idx = species * N_bins + idx;
			fvec[matrix_idx] = rho[matrix_idx] - rho_curr(idx, species) + (F_half(idx, species) - F_half(next_idx, species)) * dt / dx;
		}
	}
}

template<int dims>
BailoFiniteVolume<dims>::BailoFiniteVolume(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config) : 
		Integrator<dims>(sim_state, model, config) {
    _N_per_dim_minus_one = this->_N_per_dim - 1;
	_log2_N_per_dim = (int) std::log2(this->_N_per_dim);

	N_species = model->N_species();
    N_bins = this->_N_bins;
    fe_model = model;
    bailo<dims> = this;
    k_laplacian = this->_k_laplacian;
	dx = this->_dx;
	dt = this->_dt;
	base_csi = MultiField<double>(N_bins, fe_model->N_species());
}

template<int dims>
BailoFiniteVolume<dims>::~BailoFiniteVolume() {

}

template<int dims>
void BailoFiniteVolume<dims>::evolve() {
    static MultiField<double> solved(this->_rho.bins(), this->_model->N_species());
    rho_curr = this->_rho;

	// we compute the explicit parts (expansive contributions) at the beginning
	for(int idx = 0; idx < N_bins; idx++) {
        for(int species = 0; species < fe_model->N_species(); species++) {
    		double F_der_exp = fe_model->der_bulk_free_energy_expansive(species, rho_curr.species_view(idx));
			double interf_exp = 2 * k_laplacian * bailo<dims>->cell_laplacian(rho_curr, species, idx);
			base_csi(idx, species) = F_der_exp - 0.5 * interf_exp;
        }
    }

    int n = this->_rho.bins() * this->_model->N_species();
    int lwa = (n * (3 * n + 13)) / 2;
    static std::vector<double> wa(lwa);
    int info = fsolve(to_solve<dims>, n, this->_rho.data(), solved.data(), 1.49012e-08, wa.data(), lwa);
}

template<int dims>
void BailoFiniteVolume<dims>::_fill_coords(int coords[dims], int idx) {
	for(int d = 0; d < dims; d++) {
		coords[d] = idx & _N_per_dim_minus_one;
		idx >>= _log2_N_per_dim; // divide by N
	}
}

template<int dims>
int BailoFiniteVolume<dims>::_cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= _log2_N_per_dim; // multiply by N
	}
	return idx;
}

template<>
double BailoFiniteVolume<1>::cell_laplacian(MultiField<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + this->_N_bins) & _N_per_dim_minus_one;
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(this->_dx);
}

template<>
double BailoFiniteVolume<2>::cell_laplacian(MultiField<double> &field, int species, int idx) {
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

template class BailoFiniteVolume<1>;
template class BailoFiniteVolume<2>;

} /* namespace */
