/*
 * CahnHilliard.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CAHNHILLIARD_H_
#define SRC_CAHNHILLIARD_H_

#include "defs.h"

#include <vector>
#include <string>

namespace ch {

template<int dims>
class CahnHilliard {
public:
	int N;
	int N_minus_one;
	int bits;
	int size;
	double dt;
	double k_laplacian;
	double M;
	double H;
	FreeEnergyModel *model;

	std::vector<std::vector<double>> rho;

	CahnHilliard(FreeEnergyModel *m, cxxopts::ParseResult &options);
	~CahnHilliard();

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	double cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx);

	void evolve();
	double total_mass();

	void print_state(int species, std::ofstream &output);
	void print_density(std::string filename);
};

} /* namespace ch */

#endif /* SRC_CAHNHILLIARD_H_ */
