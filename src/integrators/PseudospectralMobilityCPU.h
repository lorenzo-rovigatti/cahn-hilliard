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

    // spectral buffers
    std::vector<std::complex<double>> rho_hat, rho_hat_copy;
    std::vector<std::complex<double>> f_der_hat;
    std::vector<std::complex<double>> mu_hat;
    std::vector<std::complex<double>> divJ_hat;

    bool use_dealias = false;
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

    // --- Option C: GMRES + tmp FFT plans/workspaces ---
    MultiField<double> tmp_real;                       // bins x species (workspace)
    std::vector<std::complex<double>> tmp_hat;         // hat_vector_size (workspace)
    std::vector<std::complex<double>> q_hat;           // hat_vector_size (workspace)
    std::vector<std::complex<double>> div_hat;         // hat_vector_size (workspace)

    fftw_plan tmp_r2c_plan = nullptr;                  // tmp_real -> tmp_hat
    fftw_plan tmp_c2r_plan = nullptr;                  // tmp_hat -> tmp_real (overwrites tmp_hat)

    // per-species reference mobility
    std::vector<double> M0_species;

    // GMRES params
    int gmres_restart = 30;
    int gmres_max_iter = 200;
    double gmres_tol = 1e-10;

    // helpers
    void _compute_M0_species();
    void _apply_A(const std::vector<double> &x, std::vector<double> &Ax);
    void _apply_Pinv(const std::vector<double> &r, std::vector<double> &z);
    int  _gmres_solve(const std::vector<double> &b, std::vector<double> &x);
    double _dot(const std::vector<double>& a, const std::vector<double>& b) const;
    double _norm2(const std::vector<double>& a) const;

    bool _use_gmres = true;

    /**
     * @brief Semi-implicit pseudospectral time step with mobility splitting.
     *
     * This method advances the Cahn–Hilliard equation with spatially varying mobility
     *
     *     ∂_t ρ = ∇·( M(x) ∇μ ),
     *     μ = μ_bulk(ρ) − 2k ∇²ρ ,
     *
     * using a *mobility-splitting* strategy:
     *
     *     M(x) = M0 + (M(x) − M0),
     *
     * where M0 is a constant reference mobility (typically the spatial maximum).
     *
     * The time discretization is Eyre-stabilized and semi-implicit:
     *
     *  • The stiff constant-mobility part
     *
     *        M0 ∇²( μ_bulk(ρ) − Sρ − 2k ∇²ρ )
     *
     *    is treated implicitly in Fourier space, yielding a diagonal k-space update
     *    with unconditional linear stability for sufficiently large S.
     *
     *  • The variable-mobility remainder
     *
     *        ∇·( (M − M0) ∇μ )
     *
     *    is treated explicitly via real-space fluxes and pseudospectral derivatives.
     *
     * In practice, this leads to the update
     *
     *     ρ^{n+1} = [ ρ^n
     *                + Δt ∇·( (M − M0) ∇μ^n )
     *                − Δt M0 ∇²( μ_bulk(ρ^n) − Sρ^n ) ]
     *               / [ 1 + Δt M0 ( S k² + 2k k⁴ ) ],
     *
     * where all Fourier-space divisions are diagonal and inexpensive.
     *
     * Properties:
     *  • Efficient: one FFT-based update per timestep, no linear solves.
     *  • Reduces to the standard semi-implicit pseudospectral scheme when M is constant.
     *  • Allows larger timesteps than fully explicit schemes.
     *
     * Limitations:
     *  • The variable-mobility interfacial contribution
     *
     *        ∇·( (M − M0) ∇( −2k ∇²ρ ) )
     *
     *    is treated explicitly; for strongly varying mobility or large k, this may
     *    affect energy-dissipation accuracy and lead to energy traces that do not
     *    closely shadow fully implicit or finite-volume schemes unless Δt is small
     *    and dealiasing is carefully applied.
     *
     * This method is best viewed as a fast, approximate solver suitable when
     * mobility variations are mild or when qualitative dynamics are sufficient.
     */
    void _evolve_simple();

    /**
     * @brief Fully implicit pseudospectral time step with variable mobility.
     *
     * This method advances the Cahn–Hilliard equation with spatially varying mobility
     *
     *     ∂_t ρ = ∇·( M(x) ∇μ ),
     *     μ = μ_bulk(ρ) − 2k ∇²ρ ,
     *
     * using a *fully implicit*, Eyre-stabilized discretization in which the variable
     * mobility multiplies the entire interfacial operator.
     *
     * The time discretization reads:
     *
     *     (ρ^{n+1} − ρ^n)/Δt
     *       = ∇·( M(x) ∇[ μ_bulk(ρ^n) − Sρ^n
     *                     + Sρ^{n+1} − 2k ∇²ρ^{n+1} ] ),
     *
     * which leads to the linear system
     *
     *     A(ρ^{n+1}) = b,
     *
     * with
     *
     *     A(x) = x − Δt ∇·( M(x) ∇( Sx − 2k ∇²x ) ),
     *     b    = ρ^n + Δt ∇·( M(x) ∇( μ_bulk(ρ^n) − Sρ^n ) ).
     *
     * The operator A is a variable-coefficient fourth-order elliptic operator and is
     * not diagonal in Fourier space. It is applied in a *matrix-free* manner using
     * pseudospectral derivatives and real-space multiplications.
     *
     * The linear system is solved at each timestep using restarted GMRES with a
     * physics-based FFT-diagonal preconditioner:
     *
     *     P(x) = x − Δt M0 ∇²( Sx − 2k ∇²x ),
     *
     * where M0 is a constant reference mobility (typically the spatial maximum).
     * The preconditioner inversion is performed exactly in Fourier space.
     *
     * Properties:
     *  • Correctly treats the full variable-mobility operator, including interfacial
     *    contributions, without approximation.
     *  • Produces energy-dissipation behavior that closely shadows finite-volume,
     *    finite-element, and fully implicit reference schemes.
     *  • Unconditionally stable with respect to the linear stiff terms.
     *
     * Costs and considerations:
     *  • Requires multiple FFTs per GMRES iteration; computationally more expensive
     *    than _evolve_simple().
     *  • Mobility is treated in a lagged fashion M = M(ρ^n); for mobility depending
     *    strongly on ρ, an outer Picard or Newton iteration may be required.
     *  • Krylov tolerance controls mass conservation and energy accuracy.
     *
     * This method should be used when quantitative accuracy of the dynamics and
     * energy dissipation is required, especially in the presence of strong spatial
     * variations in mobility.
     */
    void _evolve_full();
};

} /* namespace ch */

#endif /* PSEUDOSPECTRALMOBILITYCPU_H */
