#ifndef EULERMOBILITYCPU_H
#define EULERMOBILITYCPU_H

#include "EulerCPU.h"

#include <random>

namespace ch {

template<int dims>
class EulerMobilityCPU : public EulerCPU<dims> {
public:
    EulerMobilityCPU(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config);

    ~EulerMobilityCPU();

    void evolve() override;

    GET_NAME(EulerMobilityCPU)

protected:
    bool _supports_nonconstant_mobility() const override {
        return true;
    }

private:
    bool _with_noise = false;
    double _noise_factor = 0.0;
    std::mt19937 _generator;
    
};

} /* namespace ch */

#endif /* EULERMOBILITYCPU_H */
