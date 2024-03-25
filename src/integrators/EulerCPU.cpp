#include "EulerCPU.h"

namespace ch {

template<int dims>
EulerCPU<dims>::EulerCPU(FreeEnergyModel *model, toml::table &config) : Integrator<dims>(model, config) {
    _N_bins_minus_one = this->_N_bins - 1;
	_log2_N_bins = (int) std::log2(this->_N_bins);
}

template<int dims>
EulerCPU<dims>::~EulerCPU() {

}

template<int dims>
void EulerCPU<dims>::evolve() {
    static RhoMatrix<double> rho_der(this->_rho.bins(), this->_model->N_species());
    // we first evaluate the time derivative for all the fields
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_model->N_species(); species++) {
            rho_der(idx, species) = this->_model->der_bulk_free_energy(species, this->_rho.rho_species(idx)) - 2 * this->_k_laplacian * _cell_laplacian(this->_rho, species, idx);
        }
    }

    // and then we integrate them
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_model->N_species(); species++) {
            this->_rho(idx, species) += this->_M * _cell_laplacian(rho_der, species, idx) * this->_dt;
        }
    }
}

template<int dims>
void EulerCPU<dims>::_fill_coords(int coords[dims], int idx) {
	for(int d = 0; d < dims; d++) {
		coords[d] = idx & _N_bins_minus_one;
		idx >>= _log2_N_bins; // divide by N
	}
}

template<int dims>
int EulerCPU<dims>::_cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= _log2_N_bins; // multiply by N
	}
	return idx;
}

template<>
double EulerCPU<1>::_cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + this->_N_bins) & _N_bins_minus_one;
	int idx_p = (idx + 1) & _N_bins_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(this->_dx);
}

template<>
double EulerCPU<2>::_cell_laplacian(RhoMatrix<double> &field, int species, int idx) {
	int coords_xy[2];
	_fill_coords(coords_xy, idx);

	int coords_xmy[2] = {
			(coords_xy[0] - 1 + this->_N_bins) & _N_bins_minus_one,
			coords_xy[1]
	};

	int coords_xym[2] = {
			coords_xy[0],
			(coords_xy[1] - 1 + this->_N_bins) & _N_bins_minus_one
	};

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & _N_bins_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & _N_bins_minus_one
	};

	return (
			field(_cell_idx(coords_xmy), species) +
			field(_cell_idx(coords_xpy), species) +
			field(_cell_idx(coords_xym), species) +
			field(_cell_idx(coords_xyp), species) -
			4 * field(idx, species))
			/ SQR(this->_dx);
}

template class EulerCPU<1>;
template class EulerCPU<2>;

} /* namespace ch */
