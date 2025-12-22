#include "EulerMobilityCPU.h"

#include "../utils/utility_functions.h"

namespace ch {

template<int dims>
EulerMobilityCPU<dims>::EulerMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
		EulerCPU<dims>(sim_state, model, config) {
	std::string mobility = this->template _config_value<std::string>(config, "mobility.type");
	_with_noise = this->template _config_optional_value<bool>(config, "mobility.with_noise", false);
	
	if(_with_noise) {
		double noise_rescale_factor = this->template _config_optional_value<double>(config, "mobility.noise_rescale_factor", 1.0);
		_noise_factor = std::sqrt(2.0 / (pow_dims<dims>(this->_dx) * this->_dt)) * noise_rescale_factor;
		this->info("Integrating the Cahn-Hilliard equation with non-constant mobility and noise (noise_factor = {})", _noise_factor);

		long long int seed = this->template _config_optional_value<long long int>(config, "seed", std::time(NULL));
		_generator.seed(seed);
	}
	else {
		this->info("Integrating the Cahn-Hilliard equation with non-constant mobility and without noise");
	}
}

template<int dims>
EulerMobilityCPU<dims>::~EulerMobilityCPU() {

}

template<int dims>
void EulerMobilityCPU<dims>::evolve() {
    static MultiField<double> rho_der(this->_rho.bins(), this->_model->N_species());
	static MultiField<Gradient<dims>> flux(this->_rho.bins(), this->_model->N_species());
	static MultiField<Gradient<dims>> stochastic_flux(this->_rho.bins(), this->_model->N_species());
	static std::normal_distribution<double> normal_dist(0.0, 1.0);

    // we first evaluate the time derivative for all the fields
    this->_model->der_bulk_free_energy(this->_rho, rho_der);
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
        for(int species = 0; species < this->_model->N_species(); species++) {
            rho_der(idx, species) -= 2 * this->_k_laplacian * this->_cell_laplacian(this->_rho, species, idx);
        }
    }

	// here we use a staggered grid discretisation to avoid numerical artifacts when computing the gradients
	for(unsigned int idx = 0; idx < this->_N_bins - 1; idx++) {
		int idx_p = (idx + 1) & this->_N_per_dim_minus_one;
        for(int species = 0; species < this->_model->N_species(); species++) {
			double M_idx = this->_sim_state.mobility(idx, species);
			double M_p = this->_sim_state.mobility(idx_p, species);
			// double M_idx = this->_M * this->_rho(idx, species) / (this->_rho(idx, species) + _rho_min);
			// double M_p = this->_M * this->_rho(idx_p, species) / (this->_rho(idx_p, species) + _rho_min);
			double M_flux = 0.5 * (M_idx + M_p);
			flux(idx, species) = M_flux * this->_cell_gradient(rho_der, species, idx);

			if(_with_noise) {
				double noise_amplitude = std::sqrt(M_flux) * _noise_factor;
				for (std::size_t d = 0; d < dims; ++d) {
					stochastic_flux(idx, species)[d] = noise_amplitude * normal_dist(_generator);
				}
			}
        }
    }

	for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
		int idx_m = (idx - 1 + this->_N_bins) & this->_N_per_dim_minus_one;
        for(int species = 0; species < this->_model->N_species(); species++) {
			this->_rho(idx, species) += (this->_divergence(flux, species, idx) + this->_divergence(stochastic_flux, species, idx)) * this->_dt;
        }
    }
}

template class EulerMobilityCPU<1>;
template class EulerMobilityCPU<2>;

} /* namespace ch */
