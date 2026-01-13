/*
 * SimpleWertheim.h
 *
 *  Created on: Jul 22, 2023
 *      Author: lorenzo
 */

#ifndef SRC_MODELS_SIMPLEWERTHEIM_H_
#define SRC_MODELS_SIMPLEWERTHEIM_H_

#include "FreeEnergyModel.h"

namespace ch {

class SimpleWertheim final: public FreeEnergyModel {
public:
	SimpleWertheim(toml::table &config);
	virtual ~SimpleWertheim();
	SimpleWertheim(const SimpleWertheim &other) = default;
	SimpleWertheim(SimpleWertheim &&other) = default;

	int N_species() override {
		return 1;
	}

	void der_bulk_free_energy(field_type *rho, float *rho_der, int grid_size) override;
	void der_bulk_free_energy(const MultiField<double> &rho, MultiField<double> &rho_der) override;
	double der_bulk_free_energy_expansive(int species, const SpeciesView<double> &) override;
	double der_bulk_free_energy_contractive(int species, const SpeciesView<double> &) override;
	double bulk_free_energy(const SpeciesView<double> &) override;

	void set_mobility(field_type *rho, double M0, field_type *mobility, int grid_size) override;
	void set_mobility(const MultiField<double> &rho, double M0, MultiField<double> &mobility) override;

	GET_NAME(Simple Wertheim free energy model)

private:
	double _B2;
	int _valence;
	double _delta, _two_valence_delta;
	double _regularisation_delta, _log_delta;

	double _X(double rho);

	inline double _regularised_log(double rho) {
		return (rho < _regularisation_delta) ? _log_delta + (rho - _regularisation_delta) / (2.0 * _regularisation_delta) : std::log(rho);
	}
};

} /* namespace ch */

#endif /* SRC_MODELS_SIMPLEWERTHEIM_H_ */
