/*
 * CahnHilliard.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CAHNHILLIARD_H_
#define SRC_CAHNHILLIARD_H_

#include "defs.h"
#include "SimulationState.h"
#include "utils/Gradient.h"
#include "utils/MultiField.h"
#include "models/FreeEnergyModel.h"
#include "integrators/Integrator.h"
#include "mobility/IMobility.h"

#include <vector>
#include <string>

namespace ch {

template<int dims>
class CahnHilliard : public Object {
public:
	int N = 0;
	int N_minus_one = 0;
	int bits = 0;
	int grid_size = 0;
	double dt = 0.0;
	std::vector<double> k_laplacian;
	double dx = 0.0;
	double V_bin;
	FreeEnergyModel *model = nullptr;
	std::unique_ptr<ch::Integrator<dims>> integrator;
	std::unique_ptr<IMobility<dims>> mobility = nullptr;

	CahnHilliard(SimulationState<dims> &sim_state, FreeEnergyModel *m, toml::table &config);
	~CahnHilliard();

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	Gradient<dims> gradient(MultiField<double> &field, int species, int idx);

	void evolve();

	double average_mass();
	double average_free_energy();
	double average_pressure();
	MultiField<double> pressure();

	GET_NAME(Simulation manager)

private:
	SimulationState<dims> &_sim_state;
	double _user_to_internal, _internal_to_user;
	bool _output_ready = false;
	int _d_vec_size;
	std::string _grid_size_str;
};

} /* namespace ch */

#endif /* SRC_CAHNHILLIARD_H_ */
