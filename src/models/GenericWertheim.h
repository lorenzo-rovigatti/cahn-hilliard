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

using DeltaMap = std::map<std::pair<int, int>, double>;

struct Species {
	int idx;
	int N_patches;
	std::vector<int> patches;
};

class GenericWertheim final: public FreeEnergyModel {
public:
	GenericWertheim(toml::table &config);
	virtual ~GenericWertheim();
	GenericWertheim(const GenericWertheim &other) = default;
	GenericWertheim(GenericWertheim &&other) = default;

	int N_species() override {
		return _species.size();
	}

	double bonding_free_energy(const std::vector<double> &);
	double bulk_free_energy(const std::vector<double> &) override;

	void der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) override;
	void der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) override;
	double der_bulk_free_energy(int species, const std::vector<double> &) override;

	GET_NAME("Generic Wertheim free energy")

private:
	std::vector<Species> _species;
	DeltaMap _delta;
	int _N_patches = 0;
	double _B2, _B3 = 0;

	void _update_X(const std::vector<double> &, std::vector<double> &);
	double _der_contribution(const std::vector<double> &, int);
};

} /* namespace ch */

#endif /* SRC_MODELS_GENERICWERTHEIM_H_ */
