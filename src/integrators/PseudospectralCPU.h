/*
 * PseudospectralCPU.h
 *
 * Created on: 3/27/2024
 *     Author: Lorenzo
*/

#ifndef PSEUDOSPECTRALCPU_H
#define PSEUDOSPECTRALCPU_H

#include "Integrator.h"

#include <array>
#include <complex>
#include <fftw3.h>

namespace ch {

template<int dims>
class PseudospectralCPU : public Integrator<dims> {
public:
    PseudospectralCPU(FreeEnergyModel *model, toml::table &config);

    ~PseudospectralCPU();

    void set_initial_rho(RhoMatrix<double> &r) override;

    void evolve() override;

    GET_NAME(PseudospectralCPU)

private:
    std::array<int, dims> _reciprocal_n; // the dimensions of the grids to be transformed
	int hat_grid_size; // n1 x n2 x ... x (n_d / 2 + 1)
	int hat_vector_size; // hat_grid_size * N_species
	std::vector<std::complex<double>> rho_hat, rho_hat_copy, f_der_hat;
	RhoMatrix<double> f_der;
	std::vector<double> sqr_wave_vectors, dealiaser;

	fftw_plan rho_inverse_plan, f_der_plan;
};

} /* namespace ch */

#endif /* PSEUDOSPECTRALCPU_H */
