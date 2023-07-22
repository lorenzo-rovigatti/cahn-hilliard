/*
 * CahnHilliard.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CAHNHILLIARD_H_
#define SRC_CAHNHILLIARD_H_

#include "defs.h"
#include "models/FreeEnergyModel.h"

#include <vector>
#include <string>

namespace ch {

template<int dims>
class CahnHilliard : public Object {
public:
	int N = 0;
	int N_minus_one = 0;
	int bits = 0;
	int size = 0;
	double dt = 0.0;
	double k_laplacian = 0.0;
	double M = 0.0;
	double H = 0.0;
	FreeEnergyModel *model = nullptr;

	std::vector<std::vector<double>> rho;

	CahnHilliard(FreeEnergyModel *m, toml::table &config);
	~CahnHilliard();

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	double cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx);

	void evolve();
	double total_mass();

	void print_state(int species, std::ofstream &output);
	void print_density(std::string filename);

	GET_NAME(Simulation manager)

private:
};

} /* namespace ch */

#endif /* SRC_CAHNHILLIARD_H_ */
