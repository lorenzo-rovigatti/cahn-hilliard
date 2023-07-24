/*
 * SalehWertheim.h
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_SALEHWERTHEIM_H_
#define SRC_MODELS_SALEHWERTHEIM_H_

#include "FreeEnergyModel.h"

namespace ch {

class SalehWertheim final: public FreeEnergyModel {
public:
	SalehWertheim(toml::table &config);
	virtual ~SalehWertheim();
	SalehWertheim(const SalehWertheim &other) = default;
	SalehWertheim(SalehWertheim &&other) = default;

	int N_species() override {
		return 3;
	}

	double bonding_free_energy(std::vector<double> &);
	double bulk_free_energy(std::vector<double> &) override;
	double der_bulk_free_energy(int species, std::vector<double> &) override;

	GET_NAME("Saleh's system Wertheim free energy")

private:
	std::vector<int> _valence;
	double _B2;
	double _delta_AA, _delta_BB;
};

} /* namespace ch */

#endif /* SRC_MODELS_SALEHWERTHEIM_H_ */
