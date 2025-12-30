/*
 * PseudospectralMobilityCPU.h
 *
 * Created on: 12/30/2025
 *      Author: Lorenzo
*/

#ifndef PSEUDOSPECTRALMOBILITYCPU_H
#define PSEUDOSPECTRALMOBILITYCPU_H

#include "Integrator.h"

#include <complex>
#include <fftw3.h>

namespace ch {

/*
 * PseudospectralMobilityCPU.cpp
 *
 * Variable-mobility pseudospectral integrator:
 *   d rho / dt = div( M(x) grad mu )
 *   mu = d f_bulk / d rho - kappa laplacian rho
 *
 * Keeps k^4 term semi-implicit using a constant reference mobility M0 (max over space),
 * while computing variable-mobility flux divergence explicitly in the numerator.
 */
template<int dims>
class PseudospectralMobilityCPU : public Integrator<dims> {
public:
    PseudospectralMobilityCPU(
        SimulationState &sim_state,
        FreeEnergyModel *model,
        toml::table &config
    );

    ~PseudospectralMobilityCPU();

    void evolve() override;

protected:
    bool _supports_nonconstant_mobility() const override { return true; }

private:
    // spectral sizes
    std::array<int, dims> _reciprocal_n{};
    int _hat_grid_size = 0;
    int _hat_vector_size = 0;

    // stabilization parameter
    double _S = 0.0;
    double _M0 = 0.0;
    bool use_dealias = false;

    // spectral buffers
    std::vector<std::complex<double>> rho_hat, rho_hat_copy;
    std::vector<std::complex<double>> f_der_hat;
    std::vector<std::complex<double>> mu_hat;
    std::vector<std::complex<double>> divJ_hat;

    std::vector<double> sqr_wave_vectors, dealiaser;
    std::array<std::vector<double>, dims> kcomp;

    // real-space fields
    MultiField<double> f_der;
    std::array<MultiField<double>, dims> grad_mu;
    std::array<MultiField<double>, dims> flux;

    // spectral flux
    std::array<std::vector<std::complex<double>>, dims> flux_hat;
    std::array<std::vector<std::complex<double>>, dims> grad_mu_hat;

    // FFTW plans
    fftw_plan f_der_plan = nullptr;
    fftw_plan rho_inverse_plan = nullptr;
    std::array<fftw_plan, dims> grad_mu_inverse_plan{};
    std::array<fftw_plan, dims> flux_plan{};
};

} /* namespace ch */

#endif /* PSEUDOSPECTRALMOBILITYCPU_H */
