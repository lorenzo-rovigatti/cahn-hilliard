/*
 * GenericWertheim.h
 *
 *  Created on: Feb 19, 2025
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_GENERICWERTHEIM_H_
#define SRC_MODELS_GENERICWERTHEIM_H_

#include "FreeEnergyModel.h"

#include <set>

namespace ch {

struct Species {
	int idx;
	int N_unique_patches;
	// all the patches (e.g. [0, 0, 1])
	std::vector<int> patches;
	// only unique patches ([0, 1] for the example above)
	std::vector<int> unique_patches;
	// the multiplicity of each unique patch ([2, 1] for the example above)
	std::vector<int> unique_patch_multiplicity;
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
	std::vector<int> _unique_patches;
	// list of species with which each unique patch can interact
	std::vector<std::vector<int>> _interacting_species;
	std::vector<double> _delta;
	int _N_patches = 0;
	double _B2, _B3 = 0;

	void _update_X(const std::vector<double> &, std::vector<double> &);
	double _der_contribution(const std::vector<double> &, int);
};

} /* namespace ch */

#endif /* SRC_MODELS_GENERICWERTHEIM_H_ */
