/*
 * defs.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#define SQR(X) ((X) * (X))

#define TOML_EXCEPTIONS 0
//#define TOML_ENABLE_FORMATTERS 0
#include <toml++/toml.hpp>

namespace ch {

class FreeEnergyModel {
public:
	FreeEnergyModel(toml::table &config) {

	}

	virtual ~FreeEnergyModel() {

	}

	virtual int N_species() = 0;
	virtual double der_bulk_free_energy(int species, std::vector<double> &) = 0;
	virtual double bulk_free_energy(int species, std::vector<double> &) = 0;
};

}

#endif /* SRC_DEFS_H_ */
