#ifndef EULERMOBILITYCPU_H
#define EULERMOBILITYCPU_H

#include "EulerCPU.h"

namespace ch {

template<int dims>
class EulerMobilityCPU : public EulerCPU<dims> {
public:
    EulerMobilityCPU(FreeEnergyModel *model, toml::table &config);

    ~EulerMobilityCPU();

    void evolve() override;

    GET_NAME(EulerMobilityCPU)

protected:
    std::array<double, dims> _cell_gradient(RhoMatrix<double> &field, int species, int idx);

private:
    double _rho_min;
};

} /* namespace ch */

#endif /* EULERMOBILITYCPU_H */
