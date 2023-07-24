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

	}

	virtual ~FreeEnergyModel() {

	}

	virtual int N_species() = 0;
	virtual double der_bulk_free_energy(int species, std::vector<double> &) = 0;
	virtual double bulk_free_energy(std::vector<double> &) = 0;
};

} /* namespace ch */

#endif /* SRC_MODELS_FREEENERGYMODEL_H_ */
