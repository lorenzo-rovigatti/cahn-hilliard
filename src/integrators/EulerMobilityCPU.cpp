#include "EulerMobilityCPU.h"

#include "../utils/utility_functions.h"

namespace ch {

template<int dims>
EulerMobilityCPU<dims>::EulerMobilityCPU(FreeEnergyModel *model, toml::table &config) : EulerCPU<dims>(model, config) {
	std::string mobility = this->template _config_value<std::string>(config, "mobility.type");
	if(mobility != "regularised") {
		this->critical("The only supported non-constant mobility is 'regularised'");
	}

	_rho_min = this->template _config_value<double>(config, "mobility.rho_min");
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
    static RhoMatrix<double> rho_der(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<Gradient<dims>> flux(this->_rho.bins(), this->_model->N_species());
	static RhoMatrix<Gradient<dims>> stochastic_flux(this->_rho.bins(), this->_model->N_species());
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
			double M_idx = this->_M * this->_rho(idx, species) / (this->_rho(idx, species) + _rho_min);
			double M_p = this->_M * this->_rho(idx_p, species) / (this->_rho(idx_p, species) + _rho_min);
			double M_flux = 0.5 * (M_idx + M_p);
			flux(idx, species) = M_flux * _cell_gradient(rho_der, species, idx);

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
			this->_rho(idx, species) += (_divergence(flux, species, idx) + _divergence(stochastic_flux, species, idx)) * this->_dt;
			// printf("%lf %lf\n", _divergence(flux, species, idx), _divergence(stochastic_flux, species, idx));
        }
    }
}

template<>
Gradient<1> EulerMobilityCPU<1>::_cell_gradient(RhoMatrix<double> &field, int species, int idx) {
	int idx_p = (idx + 1) & _N_per_dim_minus_one;

	return Gradient<1>({(field(idx_p, species) - field(idx, species)) / _dx});
}

template<>
Gradient<2> EulerMobilityCPU<2>::_cell_gradient(RhoMatrix<double> &field, int species, int idx) {
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
double EulerMobilityCPU<1>::_divergence(RhoMatrix<Gradient<1>> &flux, int species, int idx) {
	int idx_m = (idx - 1 + _N_bins) & _N_per_dim_minus_one;
	return (flux(idx, species)[0] - flux(idx_m, species)[0]) / this->_dx;
}

template<int dims>
double EulerMobilityCPU<dims>::_divergence(RhoMatrix<Gradient<dims>> &flux, int species, int idx) {
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

template class EulerMobilityCPU<1>;
template class EulerMobilityCPU<2>;

} /* namespace ch */
