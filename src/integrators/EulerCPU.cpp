#include "EulerCPU.h"

namespace ch {

template<int dims>
EulerCPU<dims>::EulerCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
		Integrator<dims>(sim_state, model, config) {
    _N_per_dim_minus_one = this->_N_per_dim - 1;
	_log2_N_per_dim = (int) std::log2(this->_N_per_dim);
}

template<int dims>
EulerCPU<dims>::~EulerCPU() {

}

template<int dims>
void EulerCPU<dims>::evolve() {
    static MultiField<double> rho_der(this->_rho.bins(), this->_N_species);
    // we first evaluate the time derivative for all the fields
	this->_model->der_bulk_free_energy(this->_rho, rho_der);
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_N_species; species++) {
            rho_der(idx, species) -= 2 * this->_k_laplacian * _cell_laplacian(this->_rho, species, idx);
        }
    }

    // and then we integrate them
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_N_species; species++) {
			double total_derivative = this->_M * _cell_laplacian(rho_der, species, idx);
            this->_rho(idx, species) += total_derivative * this->_dt;
        }
    }
}

template<int dims>
void EulerCPU<dims>::_fill_coords(int coords[dims], int idx) {
	for(int d = 0; d < dims; d++) {
		coords[d] = idx & _N_per_dim_minus_one;
		idx >>= _log2_N_per_dim; // divide by N
	}
}

template<int dims>
int EulerCPU<dims>::_cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= _log2_N_per_dim; // multiply by N
	}
	return idx;
}

template<>
double EulerCPU<1>::_cell_laplacian(MultiField<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + this->_N_bins) & _N_per_dim_minus_one;
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return (field(idx_m, species) + field(idx_p, species) - 2.0 * field(idx, species)) / SQR(this->_dx);
}

template<>
double EulerCPU<2>::_cell_laplacian(MultiField<double> &field, int species, int idx) {
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

template<>
Gradient<1> EulerCPU<1>::_cell_gradient(MultiField<double> &field, int species, int idx) {
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return Gradient<1>({(field(idx_p, species) - field(idx, species)) / _dx});
}

template<>
Gradient<2> EulerCPU<2>::_cell_gradient(MultiField<double> &field, int species, int idx) {
	int coords_xy[2];
	_fill_coords(coords_xy, idx);

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & _N_per_dim_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & _N_per_dim_minus_one
	};

	return Gradient<2>({
		(field(_cell_idx(coords_xpy), species) - field(_cell_idx(coords_xy), species)) / this->_dx, 
		(field(_cell_idx(coords_xyp), species) - field(_cell_idx(coords_xy), species)) / this->_dx
	});
}

template<>
double EulerCPU<1>::_divergence(MultiField<Gradient<1>> &flux, int species, int idx) {
	int idx_m = (idx - 1 + _N_bins) & _N_per_dim_minus_one;
	return (flux(idx, species)[0] - flux(idx_m, species)[0]) / this->_dx;
}

template<int dims>
double EulerCPU<dims>::_divergence(MultiField<Gradient<dims>> &flux, int species, int idx) {
	double res = 0;
	int coords[dims], coords_m[dims];
	this->_fill_coords(coords, idx);
	memcpy(coords_m, coords, sizeof(int) * dims);

	for(int d = 0; d < dims; d++) {
		coords_m[d] = (coords[d] - 1 + this->_N_bins) & this->_N_per_dim_minus_one;
		int idx_m = this->_cell_idx(coords_m);
		res += (flux(idx, species)[d] - flux(idx_m, species)[d]) / this->_dx;
		coords_m[d] = coords[d];
	}

	return res;
}

template class EulerCPU<1>;
template class EulerCPU<2>;

} /* namespace ch */
