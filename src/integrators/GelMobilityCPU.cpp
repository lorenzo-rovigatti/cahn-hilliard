#include "GelMobilityCPU.h"

#include "../utils/utility_functions.h"

namespace ch {

template<int dims>
GelMobilityCPU<dims>::GelMobilityCPU(FreeEnergyModel *model, toml::table &config) : EulerCPU<dims>(model, config) {
	_with_noise = this->template _config_optional_value<bool>(config, "mobility.with_noise", false);
	
	if(_with_noise) {
		double noise_rescale_factor = this->template _config_optional_value<double>(config, "mobility.noise_rescale_factor", 1.0);
		_noise_factor = std::sqrt(2.0 / (pow_dims<dims>(this->_dx) * this->_dt)) * noise_rescale_factor;
		this->info("Integrating the Cahn-Hilliard equation with gel-like mobility and noise (noise_factor = {})", _noise_factor);

		long long int seed = this->template _config_optional_value<long long int>(config, "seed", std::time(NULL));
		_generator.seed(seed);
	}
	else {
		this->info("Integrating the Cahn-Hilliard equation with gel-like mobility and without noise");
	}

	double epsilon = this->template _config_value<double>(config, "landau.epsilon");
	double beta_delta_F = 10.0 / (1 - epsilon) - std::log(24000.0);
	_p_gel = exp(beta_delta_F) / (1.0 + exp(beta_delta_F));

	this->info("p_gel = {}, phi_critical = {}, c_0 = {}, M_c = {}", _p_gel, _phi_critical, _c_0, _M_c);
}

template<int dims>
GelMobilityCPU<dims>::~GelMobilityCPU() {

}

template<int dims>
void GelMobilityCPU<dims>::set_initial_rho(RhoMatrix<double> &r) {
	EulerCPU<dims>::set_initial_rho(r);
	_gel_OP = RhoMatrix<double>(this->_rho.bins(), 1);

	std::uniform_real_distribution<double> dist(0.0, 1.0);
	for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
		_gel_OP(idx, 0) = dist(_generator) * 1e-4;
	}
}

template<int dims>
void GelMobilityCPU<dims>::evolve() {
    static RhoMatrix<double> rho_der(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<Gradient<dims>> flux(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<Gradient<dims>> stochastic_flux(this->_rho.bins(), this->_model->N_species());
	static std::normal_distribution<double> normal_dist(0.0, 1.0);

    // we first evaluate the time derivative for all the fields
    this->_model->der_bulk_free_energy(this->_rho, rho_der);
    for(unsigned int idx = 0; idx < this->_N_bins; idx++) {
		double rho_tot = 0.;
        for(int species = 0; species < this->_model->N_species(); species++) {
            rho_der(idx, species) -= 2 * this->_k_laplacian * this->_cell_laplacian(this->_rho, species, idx);
			rho_tot += this->_rho(idx, species);
        }

		// update the gel OP
		double c = _gel_OP(idx, 0);
		double phi = (rho_tot + 1.0) / 2.0;
		double g = (_p_gel * phi - _phi_critical) / (1.0 - _phi_critical);
		double c_der = _M_c * (g * c - c * c);

		// if(idx == 0) {
		// 	printf("%lf %lf %lf %lf %lf %lf\n", c, phi, g, c_der, _p_gel, (1.0 - _phi_critical));
		// }

		_gel_OP(idx, 0) += c_der * this->_dt;
		// if(_gel_OP(idx, 0) < 0.) {
		// 	_gel_OP(idx, 0) = 0.;
		// }
    }

	// here we use a staggered grid discretisation to avoid numerical artifacts when computing the gradients
	for(unsigned int idx = 0; idx < this->_N_bins - 1; idx++) {
		int idx_p = (idx + 1) & this->_N_per_dim_minus_one;
        for(int species = 0; species < this->_model->N_species(); species++) {
			double M_idx = std::exp(-_gel_OP(idx, 0) / _c_0);
			double M_p = std::exp(-_gel_OP(idx_p, 0) / _c_0);
			// if(idx == 0) {
			// 	printf("%lf, %lf\n", _gel_OP(idx, 0), M_idx);
			// }
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

template class GelMobilityCPU<1>;
template class GelMobilityCPU<2>;

} /* namespace ch */
