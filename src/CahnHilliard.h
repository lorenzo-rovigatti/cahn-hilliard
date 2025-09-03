/*
 * CahnHilliard.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_CAHNHILLIARD_H_
#define SRC_CAHNHILLIARD_H_

#include "defs.h"
#include "utils/RhoMatrix.h"
#include "models/FreeEnergyModel.h"
#include "integrators/Integrator.h"

#include <vector>
#include <string>
#include <array>
#include <complex>
#include <fftw3.h>

namespace ch {

template<int dims>
class CahnHilliard : public Object {
public:
	int N = 0;
	int N_minus_one = 0;
	int bits = 0;
	int grid_size = 0; // size of the grid
	double dt = 0.0;
	double k_laplacian = 0.0;
	double dx = 0.0;
	double V_bin;
	FreeEnergyModel *model = nullptr;
	Integrator<dims> *integrator = nullptr;

	CahnHilliard(FreeEnergyModel *m, toml::table &config);
	~CahnHilliard();

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	std::array<double, dims> gradient(RhoMatrix<double> &field, int species, int idx);

	void evolve();

	double average_mass();
	double average_free_energy();
	double average_pressure();
	void print_species_density(int species, const std::string &filename, long long int t);
	void print_species_density(int species, std::ofstream &output, long long int t);
	void print_total_density(const std::string &filename, long long int t);
	void print_pressure(const std::string &filename, long long int t);
	void print_pressure(std::ofstream &output, long long int t);

	GET_NAME(Simulation manager)

private:
	double _user_to_internal, _internal_to_user;
	bool _output_ready = false;
	int _d_vec_size;
	std::string _grid_size_str;

	double _density_to_user(double v);
};

} /* namespace ch */

#endif /* SRC_CAHNHILLIARD_H_ */
