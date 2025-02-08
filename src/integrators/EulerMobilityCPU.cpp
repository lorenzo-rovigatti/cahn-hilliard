#include "EulerMobilityCPU.h"

namespace ch {

template<int dims>
EulerMobilityCPU<dims>::EulerMobilityCPU(FreeEnergyModel *model, toml::table &config) : EulerCPU<dims>(model, config) {
	std::string mobility = this->template _config_value<std::string>(config, "mobility.type");
	if(mobility != "regularised") {
		this->critical("The only supported non-constant mobility is 'regularised'");
	}

	_rho_min = this->template _config_value<double>(config, "mobility.rho_min");
}

template<int dims>
EulerMobilityCPU<dims>::~EulerMobilityCPU() {

}

template<int dims>
void EulerMobilityCPU<dims>::evolve() {
    static RhoMatrix<double> rho_der(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<double> M_mu(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<double> flux(this->_rho.bins(), this->_model->N_species());

    // we first evaluate the time derivative for all the fields
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_model->N_species(); species++) {
            rho_der(idx, species) = this->_model->der_bulk_free_energy(species, this->_rho.rho_species(idx)) - 2 * this->_k_laplacian * this->_cell_laplacian(this->_rho, species, idx);
        }
    }

	// here we use a staggered grid discretisation to avoid numerical artifacts when computing the gradients
	for(unsigned int idx = 0; idx < this->_N_bins - 1; idx++) {
		int idx_p = (idx + 1) & this->_N_per_dim_minus_one;
        for(int species = 0; species < this->_model->N_species(); species++) {
			double M_idx = this->_M * this->_rho(idx, species) / (this->_rho(idx, species) + _rho_min);
			double M_p = this->_M * this->_rho(idx_p, species) / (this->_rho(idx_p, species) + _rho_min);
			double M_flux = 0.5 * (M_idx + M_p);
			flux(idx, species) = M_flux * (rho_der(idx_p, species) - rho_der(idx, species)) / this->_dx;
        }
    }

	for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
		int idx_m = (idx - 1 + this->_N_bins) & this->_N_per_dim_minus_one;
        for(int species = 0; species < this->_model->N_species(); species++) {
			double divergence = (flux(idx, species) - flux(idx_m, species)) / this->_dx;
			this->_rho(idx, species) += divergence * this->_dt;
        }
    }
}

template<>
std::array<double, 1> EulerMobilityCPU<1>::_cell_gradient(RhoMatrix<double> &field, int species, int idx) {
	int idx_m = (idx - 1 + _N_bins) & _N_per_dim_minus_one;
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return std::array<double, 1>({(field(idx_p, species) - field(idx_m, species)) / (2.0 * _dx)});
}

template<>
std::array<double, 2> EulerMobilityCPU<2>::_cell_gradient(RhoMatrix<double> &field, int species, int idx) {
	critical("The regularised mobility is not supported by 2D systems yet");
	return std::array<double, 2>({0, 0});
}

template class EulerMobilityCPU<1>;
template class EulerMobilityCPU<2>;

} /* namespace ch */
