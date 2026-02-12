#include "EulerCPU.h"

namespace ch {

template<int dims>
EulerCPU<dims>::EulerCPU(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config) : 
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
    double M = this->_sim_state.mobility(0, 0); // constant mobility
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_N_species; species++) {
			double total_derivative = M * _cell_laplacian(rho_der, species, idx);
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

template<int dims>
double EulerCPU<dims>::_cell_laplacian(MultiField<double> &field, int species, int idx) {
    if constexpr (dims == 1) {
        int idx_m = (idx - 1 + this->_N_per_dim) & this->_N_per_dim_minus_one;
        int idx_p = (idx + 1) & this->_N_per_dim_minus_one;

        return (field(idx_m, species)
              + field(idx_p, species)
              - 2.0 * field(idx, species)) / SQR(this->_dx);
    } 
	else {
        int coords[dims];
        int coords_n[dims];
        this->_fill_coords(coords, idx);

        double sum = 0.0;

        for(int d = 0; d < dims; d++) {
            // minus direction
            memcpy(coords_n, coords, sizeof(coords));
            coords_n[d] = (coords[d] - 1 + this->_N_per_dim) & this->_N_per_dim_minus_one;
            sum += field(this->_cell_idx(coords_n), species);

            // plus direction
            memcpy(coords_n, coords, sizeof(coords));
            coords_n[d] = (coords[d] + 1) & this->_N_per_dim_minus_one;
            sum += field(this->_cell_idx(coords_n), species);
        }

        return (sum - 2.0 * dims * field(idx, species)) / SQR(this->_dx);
    }
}

template<int dims>
Gradient<dims> EulerCPU<dims>::_cell_gradient(MultiField<double> &field, int species, int idx) {
    Gradient<dims> grad{};

    if constexpr (dims == 1) {
        int idx_p = (idx + 1) & this->_N_per_dim_minus_one;
        grad[0] = (field(idx_p, species) - field(idx, species)) / this->_dx;
    } 
	else {
        int coords[dims];
        int coords_p[dims];
        this->_fill_coords(coords, idx);

        for(int d = 0; d < dims; d++) {
            memcpy(coords_p, coords, sizeof(coords));
            coords_p[d] = (coords[d] + 1) & this->_N_per_dim_minus_one;
            grad[d] = (field(this->_cell_idx(coords_p), species) - field(idx, species)) / this->_dx;
        }
    }

    return grad;
}

template<int dims>
double EulerCPU<dims>::_divergence(MultiField<Gradient<dims>> &flux, int species, int idx) {
	double res = 0;
	int coords[dims], coords_m[dims];
	this->_fill_coords(coords, idx);
	memcpy(coords_m, coords, sizeof(coords));

	for(int d = 0; d < dims; d++) {
		coords_m[d] = (coords[d] - 1 + this->_N_per_dim) & this->_N_per_dim_minus_one;
		int idx_m = this->_cell_idx(coords_m);
		res += (flux(idx, species)[d] - flux(idx_m, species)[d]) / this->_dx;
		coords_m[d] = coords[d];
	}

	return res;
}

template class EulerCPU<1>;
template class EulerCPU<2>;
template class EulerCPU<3>;

} /* namespace ch */
