#ifndef EULERCPU_H
#define EULERCPU_H

#include "Integrator.h"

#include "../utils/Gradient.h"

namespace ch {

template<int dims>
class EulerCPU : public Integrator<dims> {
public:
    EulerCPU(FreeEnergyModel *model, toml::table &config);

    ~EulerCPU();

    void evolve() override;

    GET_NAME(EulerCPU)

protected:
    int _N_per_dim_minus_one;
    int _log2_N_per_dim;

    bool _couple_pressure;
    double _pressure_lambda, _pressure_target;

    void _fill_coords(int coords[dims], int idx);
    int _cell_idx(int coords[dims]);
    double _cell_laplacian(RhoMatrix<double> &field, int species, int idx);
    Gradient<dims> _cell_gradient(RhoMatrix<double> &field, int species, int idx);
    double _divergence(RhoMatrix<Gradient<dims>> &gradients, int species, int idx);
};

} /* namespace ch */

#endif /* EULERCPU_H */
