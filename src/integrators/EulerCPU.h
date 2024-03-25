#ifndef EULERCPU_H
#define EULERCPU_H

#include "Integrator.h"

namespace ch {

template<int dims>
class EulerCPU : public Integrator<dims> {
public:
    EulerCPU(FreeEnergyModel *model, toml::table &config);

    ~EulerCPU();

    void evolve() override;

private:
    int _N_bins_minus_one;
    int _log2_N_bins;

    void _fill_coords(int coords[dims], int idx);
    int _cell_idx(int coords[dims]);
    double _cell_laplacian(RhoMatrix<double> &field, int species, int idx);
};

} /* namespace ch */

#endif /* EULERCPU_H */