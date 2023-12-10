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
	double M = 0.0;
	double dx = 0.0;
	double V_bin;
	FreeEnergyModel *model = nullptr;

	RhoMatrix<double> rho;

	CahnHilliard(FreeEnergyModel *m, toml::table &config);
	~CahnHilliard();

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	std::array<double, dims> gradient(RhoMatrix<double> &field, int species, int idx);
	double cell_laplacian(RhoMatrix<double> &field, int species, int idx);

	void evolve();

	double average_mass();
	double average_free_energy();
	void print_species_density(int species, const std::string &filename);
	void print_species_density(int species, std::ofstream &output);
	void print_total_density(const std::string &filename);

	GET_NAME(Simulation manager)

private:
	double _user_to_internal, _internal_to_user;
	bool _reciprocal = false;
	bool _use_CUDA = false;
	bool _output_ready = false;
	int _d_vec_size;
	RhoMatrix<field_type> _h_rho;
	field_type *_d_rho = nullptr;
	float *_d_rho_der = nullptr;

	void _evolve_direct();
	double _density_to_user(double v);

	void _init_CUDA(toml::table &config);
	void _CPU_GPU();
	void _GPU_CPU();

	// pseudospectral stuff
	std::array<int, dims> _reciprocal_n; // the dimensions of the grids to be transformed
	int hat_grid_size; // n1 x n2 x ... x (n_d / 2 + 1)
	int hat_vector_size; // hat_grid_size * N_species
	std::vector<std::complex<double>> rho_hat, rho_hat_copy, f_der_hat;
	RhoMatrix<double> f_der;
	std::vector<double> sqr_wave_vectors, dealiaser;

	fftw_plan rho_inverse_plan, f_der_plan;

#ifndef NOCUDA
	cufftFieldComplex *_d_rho_hat = nullptr, *_d_rho_hat_copy;
	cufftComplex *_d_f_der_hat = nullptr;
	float *_d_sqr_wave_vectors = nullptr; 
	float *_d_dealiaser = nullptr;

	cufftHandle _d_rho_inverse_plan, _d_f_der_plan;
#endif

	void _evolve_reciprocal();
};

} /* namespace ch */

#endif /* SRC_CAHNHILLIARD_H_ */
