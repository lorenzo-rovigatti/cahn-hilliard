/*
 * FreeEnergyModel.h
 *
 *  Created on: Jul 21, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_FREEENERGYMODEL_H_
#define SRC_MODELS_FREEENERGYMODEL_H_

#include "../Object.h"

#include <vector>

namespace ch {

class FreeEnergyModel : public Object {
public:
	FreeEnergyModel(toml::table &config) {
		_use_CUDA = _config_optional_value<bool>(config, "use_CUDA", false);
	}

	virtual ~FreeEnergyModel() {

	}

	virtual int N_species() = 0;
	virtual void der_bulk_free_energy(number *rho, number *rho_der, int grid_size) {
		critical("this model does not support CUDA simulations");
	}
	virtual double der_bulk_free_energy(int species, std::vector<double> &) = 0;
	virtual double bulk_free_energy(std::vector<double> &) = 0;

private:
	bool _use_CUDA;
};

} /* namespace ch */

#endif /* SRC_MODELS_FREEENERGYMODEL_H_ */