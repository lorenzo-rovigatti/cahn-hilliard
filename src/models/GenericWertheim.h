/*
 * GenericWertheim.h
 *
 *  Created on: Feb 19, 2025
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_GENERICWERTHEIM_H_
#define SRC_MODELS_GENERICWERTHEIM_H_

#include "FreeEnergyModel.h"

namespace ch {

class GenericWertheim final: public FreeEnergyModel {
public:
	GenericWertheim(toml::table &config);
	virtual ~GenericWertheim();
	GenericWertheim(const GenericWertheim &other) = default;
	GenericWertheim(GenericWertheim &&other) = default;

	int N_species() override {
		return 4;
	}

	double bonding_free_energy(const std::vector<double> &);
	double bulk_free_energy(const std::vector<double> &) override;

	void der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) override;
	double der_bulk_free_energy(int species, const std::vector<double> &) override;

	GET_NAME("Generic Wertheim free energy")

private:
	int _valence, _linker_partial_valence;
	double _B2, _B3 = 0;
	double _delta_AA, _delta_BB, _delta_CC;

	double _der_contribution(const std::vector<double> &, int);
};

} /* namespace ch */

#endif /* SRC_MODELS_GENERICWERTHEIM_H_ */
