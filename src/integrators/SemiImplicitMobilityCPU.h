#ifndef SEMIIMPLICITMOBILITYCPU_H
#define SEMIIMPLICITMOBILITYCPU_H

#include "PseudospectralCPU.h"

#include <random>
#include <vector>

namespace ch {

template<int dims>
class SemiImplicitMobilityCPU : public PseudospectralCPU<dims> {
public:
    SemiImplicitMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);
    ~SemiImplicitMobilityCPU();

    void evolve() override;

    GET_NAME(SemiImplicitMobilityCPU)

protected:
    bool _supports_nonconstant_mobility() const override { return true; }

private:
    double _M0 = 1.0;
    double _rho_floor = 0.0;
    int _N_per_dim_minus_one;

    bool use_dealias;

    // real-space buffers
    MultiField<double> mu_real;
    MultiField<double> corr_real;

    // spectral buffer for correction
    std::vector<std::complex<double>> corr_hat;

    // helpers
    void _fill_coords(int coords[dims], int idx) const;
    int _cell_idx(const int coords[dims]) const;
    double _cell_laplacian(MultiField<double> &field, int species, int idx) const;
    void compute_mu();
    void compute_correction(int species);
};

} /* namespace ch */

#endif /* SEMIIMPLICITMOBILITYCPU_H */
