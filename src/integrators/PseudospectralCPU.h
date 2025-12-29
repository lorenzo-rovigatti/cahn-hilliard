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

/**
 * FFT semi-implicit Cahn–Hilliard with variable mobility via mobility splitting:
 *
 *   rho_t = div(M grad(mu)),   mu = f'(rho) - 2*k*lap(rho)
 *
 * Split mobility: M = M0 + (M - M0)
 *
 * Semi-implicit scheme:
 *   (1 + dt*2*k*M0*k^4) rho_hat^{n+1}
 *     = rho_hat^n - dt*M0*k^2*f'_hat(rho^n) + dt*corr_hat
 *
 * where corr = div((M-M0) grad(mu^n)) computed in real space (FD) and FFT'd.
 *
 * Config keys:
 *   mobility.M0                  (double, default: mobility(0,0))
 *   semi_implicit.rho_floor       (double, default: 0.0)  // clamp rho before calling free energy derivative
 *   semi_implicit.dealias         (bool,   default: false)
 */
template<int dims>
class PseudospectralCPU : public Integrator<dims> {
public:
    PseudospectralCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);

    ~PseudospectralCPU();

    void evolve() override;

    GET_NAME(PseudospectralCPU)

protected:
    std::array<int, dims> _reciprocal_n; // the dimensions of the grids to be transformed
	int hat_grid_size; // n1 x n2 x ... x (n_d / 2 + 1)
	int hat_vector_size; // hat_grid_size * N_species
	std::vector<std::complex<double>> rho_hat, rho_hat_copy, f_der_hat;
	MultiField<double> f_der;
	std::vector<double> sqr_wave_vectors, dealiaser;

	fftw_plan rho_inverse_plan, f_der_plan;
};

} /* namespace ch */

#endif /* PSEUDOSPECTRALCPU_H */
