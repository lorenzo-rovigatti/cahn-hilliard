#include "Integrator.h"

namespace ch {

template<int dims>
Integrator<dims>::Integrator(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config) : 
                _sim_state(sim_state),
                _rho(sim_state.rho),
                _model(model) {

    _mobility_type = _config_optional_value<std::string>(config, "mobility.type", "constant");
    _k_laplacian = _config_array_values<double>(config, "k", model->N_species());
    _dt = _config_value<double>(config, "dt");
    _dx = _config_optional_value<double>(config, "dx", 1.0);

    _dx *= sim_state.user_to_internal; // proportional to m
    for(auto &k : _k_laplacian) {
        k *= std::pow(sim_state.user_to_internal, 5); // proportional to m^5
    }

    _N_per_dim = _config_value<int>(config, "N");
    _N_bins = _N_per_dim;
	for(int i = 1; i < dims; i++) {
		_N_bins *= _N_per_dim;
	}
    _N_species = model->N_species();

    std::stringstream k_laplacian_str;
    for(auto k : _k_laplacian) {
        k_laplacian_str << k << " ";
    }

    info("Integrator initialized with dt = {}, dx = {}, k = {}", _dt, _dx, k_laplacian_str.str());
}

template<int dims>
Integrator<dims>::~Integrator() {

}

template<int dims>
void Integrator<dims>::validate() {
    if(!_supports_nonconstant_mobility() && _mobility_type != "constant") {
        this->critical("The selected integrator only supports constant mobility");
    }
}

template class Integrator<1>;
template class Integrator<2>;
template class Integrator<3>;

} /* namespace ch */
