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

template<int dims>
class FreeEnergyMobility: public IMobility<dims> {
public:
    FreeEnergyMobility(SimulationState<dims>& state, double M0) :
            _M0(M0 / state.user_to_internal), 
            IMobility<dims>(state) {
    }

    ~FreeEnergyMobility() {}

    void update_mobility() override {
        if(this->_sim_state.use_CUDA) {
#ifndef NOCUDA
            this->_sim_state.model->set_mobility(this->_sim_state.CUDA_rho, _M0, this->_sim_state.CUDA_mobility, this->_sim_state.rho.size());
#endif
        }
        else {
            this->_sim_state.model->set_mobility(this->_sim_state.rho, _M0, this->_sim_state.mobility);
        }
    }

    GET_NAME(FreeEnergyMobility)

private:
    double _M0;
};

} /* namespace ch */

#endif /* FREEENERGYMOBILITY_H */
