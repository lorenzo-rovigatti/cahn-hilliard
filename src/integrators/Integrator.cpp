#include "Integrator.h"

namespace ch {

template<int dims>
Integrator<dims>::Integrator(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config) : 
                _sim_state(sim_state),
                _rho(sim_state.rho),
                _model(model) {
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
    _N_species = model->N_species();

    info("Integrator initialized with dt = {}, dx = {}, k = {}, M = {}", _dt, _dx, _k_laplacian, _M);
}

template<int dims>
Integrator<dims>::~Integrator() {

}

template class Integrator<1>;
template class Integrator<2>;

} /* namespace ch */
