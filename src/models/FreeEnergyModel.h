/*
 * FreeEnergyModel.h
 *
 *  Created on: Jul 21, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_FREEENERGYMODEL_H_
#define SRC_MODELS_FREEENERGYMODEL_H_

#include "../Object.h"
#include "../utils/RhoMatrix.h"

#include <vector>

namespace ch {

class FreeEnergyModel : public Object {
public:
	FreeEnergyModel(toml::table &config) {
		_use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);
		_user_to_internal = 1.0 / _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
		_density_conversion_factor = CUB(_user_to_internal);
	}

	virtual ~FreeEnergyModel() {

	}

	virtual int N_species() = 0;
	virtual void der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size) {
		critical("this model does not support CUDA simulations");
	}
	virtual double der_bulk_free_energy_expansive(int species, const std::vector<double> &) {
		critical("this model does not support expansive/contractive splitting");
		return 0.;
	}
	virtual double der_bulk_free_energy_contractive(int species, const std::vector<double> &) {
		critical("this model does not support expansive/contractive splitting");
		return 0.;
	}
	virtual double der_bulk_free_energy(int species, const std::vector<double> &) = 0;
	virtual void der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) = 0;
	virtual double bulk_free_energy(const std::vector<double> &) = 0;

protected:
	double _user_to_internal;
	double _density_conversion_factor;
	bool _use_CUDA;
};

} /* namespace ch */

#endif /* SRC_MODELS_FREEENERGYMODEL_H_ */
