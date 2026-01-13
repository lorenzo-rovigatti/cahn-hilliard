/*
 * Landau.h
 *
 *  Created on: Jul 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_LANDAU_H_
#define SRC_MODELS_LANDAU_H_

#include "../defs.h"
#include "FreeEnergyModel.h"

namespace ch {

class Landau final: public FreeEnergyModel {
public:
	Landau(toml::table &config);
	virtual ~Landau();
	Landau(const Landau &other) = default;
	Landau(Landau &&other) = default;

	int N_species() override;

	void der_bulk_free_energy(field_type *psi, float *psi_der, int grid_size) override;
	void der_bulk_free_energy(const MultiField<double> &rho, MultiField<double> &rho_der) override;
	double bulk_free_energy(const SpeciesView<double> &) override;
	double pressure(int species, const SpeciesView<double> &) override;

	GET_NAME(Landau free energy model)

private:
	double _epsilon = 0.0;
};

} /* namespace ch */

#endif /* SRC_MODELS_LANDAU_H_ */
