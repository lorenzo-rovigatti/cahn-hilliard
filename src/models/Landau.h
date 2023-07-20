/*
 * Landau.h
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_LANDAU_H_
#define SRC_MODELS_LANDAU_H_

#include "../defs.h"

namespace ch {

class Landau: public FreeEnergyModel {
public:
	Landau(toml::table &config);
	virtual ~Landau();
	Landau(const Landau &other) = default;
	Landau(Landau &&other) = default;

	int N_species() override;
	double der_bulk_free_energy(int species, std::vector<double> &) override;
	double bulk_free_energy(int species, std::vector<double> &) override;

	GET_NAME(Landau)

private:
	double _epsilon = 0.0;
};

} /* namespace ch */

#endif /* SRC_MODELS_LANDAU_H_ */
