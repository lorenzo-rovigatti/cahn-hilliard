#ifndef GELMOBILITYCPU_H
#define GELMOBILITYCPU_H

#include "EulerCPU.h"

#include <random>

namespace ch {

template<int dims>
class GelMobilityCPU : public EulerCPU<dims> {
public:
    GelMobilityCPU(FreeEnergyModel *model, toml::table &config);

    ~GelMobilityCPU();

    void set_initial_rho(RhoMatrix<double> &r) override;

    void evolve() override;

    GET_NAME(GelMobilityCPU)

private:
    bool _with_noise = false;
    double _noise_factor = 0.0;
    std::mt19937 _generator;

    double _phi_critical = 0.5;
    double _c_0 = 0.01;
    double _M_c = 0.02;
    double _p_gel;

    RhoMatrix<double> _gel_OP;
};

} /* namespace ch */

#endif /* GELMOBILITYCPU_H */
