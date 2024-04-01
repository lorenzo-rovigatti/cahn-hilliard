/*
 * BailoFiniteVolume.cpp
 *
 * Created on: 4/1/2024
 *     Author: Lorenzo
*/

#include "BailoFiniteVolume.h"

#include <fsolve/fsolve.hpp>

namespace ch {

int N_bins;
double k_laplacian;
FreeEnergyModel *model;

template<int dims> BailoFiniteVolume<dims> *bailo;

template<int dims>
void to_solve(int n, double *x, double *fvec) {
    static RhoMatrix<double> rho_der(N_bins, model->N_species());

    memset(fvec, 0, sizeof(double) * n);

    // we first evaluate the time derivative for all the fields
    for(unsigned int idx = 0; idx < N_bins; idx++) {
        for(int species = 0; species < model->N_species(); species++) {
            int matrix_idx = species * N_bins + idx;
            // fvec[matrix_idx] += model->der_bulk_free_energy_contractive();
            rho_der(idx, species) = model->der_bulk_free_energy(species, bailo<dims>->rho().rho_species(idx)) - 2 * k_laplacian * bailo<dims>->cell_laplacian(bailo<dims>->rho(), species, idx);
        }
    }
}

template<int dims>
BailoFiniteVolume<dims>::BailoFiniteVolume(FreeEnergyModel *model, toml::table &config) : Integrator<dims>(model, config) {
    _N_per_dim_minus_one = this->_N_per_dim - 1;
	_log2_N_per_dim = (int) std::log2(this->_N_per_dim);

    N_bins = this->_N_bins;
    model = this->_model;
    bailo<dims> = this;
    k_laplacian = this->_k_laplacian;
}

template<int dims>
BailoFiniteVolume<dims>::~BailoFiniteVolume() {

}

template<int dims>
void BailoFiniteVolume<dims>::evolve() {
    static RhoMatrix<double> rho_next(this->_rho.bins(), this->_model->N_species());
    static RhoMatrix<double> solved(this->_rho.bins(), this->_model->N_species());
    rho_next = this->_rho;

    int n = this->_rho.bins() * this->_model->N_species();
    int lwa = (n*(3*n+13)) / 2;
    static std::vector<double> wa(lwa);
    fsolve(to_solve<dims>, n, this->_rho.data(), solved.data(), 1.49012e-08, wa.data(), lwa);
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
double BailoFiniteVolume<1>::cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + this->_N_bins) & _N_per_dim_minus_one;
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(this->_dx);
}

template<>
double BailoFiniteVolume<2>::cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
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
