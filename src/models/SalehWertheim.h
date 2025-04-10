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

	double bonding_free_energy(const std::vector<double> &);
	double bulk_free_energy(const std::vector<double> &) override;

	void der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) override;
	void der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) override;

	GET_NAME("Saleh's system Wertheim free energy")

private:
	std::vector<int> _valence;
	int _linker_half_valence;
	double _B2, _B3 = 0;
	double _delta_AA, _delta_BB;

	double _der_contribution(const std::vector<double> &, int);
};

} /* namespace ch */

#endif /* SRC_MODELS_SALEHWERTHEIM_H_ */
