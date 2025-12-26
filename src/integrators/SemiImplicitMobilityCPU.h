#ifndef SEMIIMPLICITMOBILITYCPU_H
#define SEMIIMPLICITMOBILITYCPU_H

#include "EulerCPU.h"

#include <random>
#include <vector>

namespace ch {

/**
 * @brief Stabilized semi-implicit (IMEX) integrator for Cahn–Hilliard with (possibly) non-constant mobility.
 *
 * Chemical potential split:
 *   mu^{n+1} = d f_bulk/d rho (rho^n) - S*rho^n + S*rho^{n+1} - 2*k_laplacian * lap(rho^{n+1}).
 *
 * Per species linear system:
 *   A rho^{n+1} = b
 * with
 *   A(x) = x - dt * [ S * D1(x) - 2*k_laplacian * D2(x) ]
 *   D1(x) = div( M^n grad(x) )
 *   D2(x) = div( M^n grad(lap(x)) )
 * and
 *   b = rho^n + dt * div( M^n grad( d f_bulk/d rho(rho^n) - S*rho^n ) )
 *       + dt * div(stochastic_flux)   (optional)
 *
 * The solver uses a matrix-free Conjugate Gradient (CG) on each species independently.
 *
 * Config keys (all optional):
 *   semi_implicit.stabilization        (double, default 0.0)
 *   semi_implicit.rho_floor            (double, default 0.0)  // clamps densities before evaluating df_bulk/drho
 *   semi_implicit.cg_tol               (double, default 1e-8)
 *   semi_implicit.cg_max_iter          (int,    default 500)
 *   mobility.with_noise                (bool,   default false)
 *   mobility.noise_rescale_factor      (double, default 1.0)
 */
template<int dims>
class SemiImplicitMobilityCPU : public EulerCPU<dims> {
public:
    SemiImplicitMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);
    ~SemiImplicitMobilityCPU();

    void evolve() override;

    GET_NAME(SemiImplicitMobilityCPU)

protected:
    bool _supports_nonconstant_mobility() const override { return true; }

private:
    // Stabilization parameter S
    double _S = 0.0;
    double _rho_floor = 0.0;

    // CG settings
    double _cg_tol = 1e-8;
    int _cg_max_iter = 500;

    // Optional stochastic flux (same discretization as EulerMobilityCPU)
    bool _with_noise = false;
    double _noise_factor = 0.0;
    std::mt19937 _generator;

    // Build div( M grad(phi) ) into out (scalar vectors size N_bins)
    void _div_M_grad_from_scalar(const std::vector<double> &phi, int species, std::vector<double> &out);

    // Discrete laplacian on scalar vector
    void _laplacian_scalar(const std::vector<double> &in, std::vector<double> &out);

    void _apply_D1(const std::vector<double> &x, int species, std::vector<double> &out);
    void _apply_D2(const std::vector<double> &x, int species, std::vector<double> &out);
    void _apply_A(const std::vector<double> &x, int species, std::vector<double> &out);

    // CG solve A x = b, returns #iters (or -1 if not converged)
    int _solve_CG(int species, const std::vector<double> &b, std::vector<double> &x);

    static double _dot(const std::vector<double> &a, const std::vector<double> &b);
};

} /* namespace ch */

#endif /* SEMIIMPLICITMOBILITYCPU_H */
