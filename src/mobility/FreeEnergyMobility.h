/*
 * FreeEnergyMobility.h
 *
 * Created on: 12/22/2025
 *      Author: Lorenzo
*/

#ifndef FREEENERGYMOBILITY_H
#define FREEENERGYMOBILITY_H

#include "IMobility.h"

namespace ch {

class FreeEnergyMobility: public IMobility {
public:
    FreeEnergyMobility(SimulationState& state, double M0) :
            _M0(M0 / state.user_to_internal), 
            IMobility(state) {
    }

    ~FreeEnergyMobility() {}

    void update_mobility() override {
        if(_sim_state.use_CUDA) {
#ifndef NOCUDA
            _sim_state.model->set_mobility(_sim_state.CUDA_rho, _M0, _sim_state.CUDA_mobility, _sim_state.rho.size());
#endif
        }
        else {
            _sim_state.model->set_mobility(_sim_state.rho, _M0, _sim_state.mobility);
        }
    }

    GET_NAME(FreeEnergyMobility)

private:
    double _M0;
};

} /* namespace ch */

#endif /* FREEENERGYMOBILITY_H */
