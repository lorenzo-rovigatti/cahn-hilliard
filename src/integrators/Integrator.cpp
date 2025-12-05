#include "Integrator.h"

namespace ch {

template<int dims>
Integrator<dims>::Integrator(FreeEnergyModel *model, toml::table &config) : _model(model) {
    _k_laplacian = _config_optional_value<double>(config, "k", 1.0);
    _dt = _config_value<double>(config, "dt");
    _M = _config_optional_value<double>(config, "M", 1.0);
    _dx = _config_optional_value<double>(config, "dx", 1.0);

    _internal_to_user = _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
    _user_to_internal = 1.0 / _internal_to_user;

    _dx *= _user_to_internal; // proportional to m
	_M /= _user_to_internal; // proportional to m^-1
	_k_laplacian *= std::pow(_user_to_internal, 5); // proportional to m^5

    _N_per_dim = _config_value<int>(config, "N");
    _N_bins = _N_per_dim;
	for(int i = 1; i < dims; i++) {
		_N_bins *= _N_per_dim;
	}
}

template<int dims>
Integrator<dims>::~Integrator() {

}

template<int dims>
void Integrator<dims>::set_initial_rho(MultiField<double> &r) {
    if(r.bins() != _N_bins) {
        critical("The size of the grid of the initial density ({}) differs from the expected value {}", r.bins(), _N_bins);
    }

    _rho = r;
}

template class Integrator<1>;
template class Integrator<2>;

} /* namespace ch */
