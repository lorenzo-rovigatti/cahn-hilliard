/*
 * BailoFiniteVolume.h
 *
 * Created on: 4/1/2024
 *      Author: Lorenzo
*/

#ifndef BAILOFINITEVOLUME_H
#define BAILOFINITEVOLUME_H

#include "Integrator.h"

namespace ch {

template<int dims>
class BailoFiniteVolume : public Integrator<dims> {
public:
    BailoFiniteVolume(SimulationState<dims> &sim_state, FreeEnergyModel *model, toml::table &config);

    ~BailoFiniteVolume();

    double cell_laplacian(MultiField<double> &field, int species, int idx);
    void evolve() override;

    GET_NAME(Bailo integrator)

private:
    int _N_per_dim_minus_one;
    int _log2_N_per_dim;

    void _fill_coords(int coords[dims], int idx);
    int _cell_idx(int coords[dims]);
};

} /* namespace */

#endif /* BAILOFINITEVOLUME_H */
