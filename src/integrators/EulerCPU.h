#ifndef EULERCPU_H
#define EULERCPU_H

#include "Integrator.h"

#include "../utils/Gradient.h"

namespace ch {

template<int dims>
class EulerCPU : public Integrator<dims> {
public:
    EulerCPU(SimulationState &sim_state, FreeEnergyModel *model, toml::table &config);

    ~EulerCPU();

    void evolve() override;

    GET_NAME(EulerCPU)

protected:
    int _N_per_dim_minus_one;
    int _log2_N_per_dim;

    void _fill_coords(int coords[dims], int idx);
    int _cell_idx(int coords[dims]);
    double _cell_laplacian(MultiField<double> &field, int species, int idx);
    Gradient<dims> _cell_gradient(MultiField<double> &field, int species, int idx);
    double _divergence(MultiField<Gradient<dims>> &gradients, int species, int idx);
};

} /* namespace ch */

#endif /* EULERCPU_H */
