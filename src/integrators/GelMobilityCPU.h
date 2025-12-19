#ifndef GELMOBILITYCPU_H
#define GELMOBILITYCPU_H

#include "EulerCPU.h"

#include <random>

namespace ch {

template<int dims>
class GelMobilityCPU : public EulerCPU<dims> {
public:
    GelMobilityCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);

    ~GelMobilityCPU();

    void evolve() override;

    GET_NAME(GelMobilityCPU)

private:
    bool _with_noise = false;
    double _noise_factor = 0.0;
    std::mt19937 _generator;

    double _phi_critical;
    double _c_0;
    double _M_c;
    double _p_gel;

    MultiField<double> _gel_OP;
};

} /* namespace ch */

#endif /* GELMOBILITYCPU_H */
