/*
 * RegularisedMobility.h
 *
 * Created on: 12/22/2025
 *      Author: Lorenzo
*/

#ifndef REGULARISEDMOBILITY_H
#define REGULARISEDMOBILITY_H

#include "IMobility.h"

namespace ch {

template<int dims>
class RegularisedMobility: public IMobility<dims> {
public:
    RegularisedMobility(SimulationState<dims>& state, double M, double rho_min) :
            _M(M / state.user_to_internal), 
            _rho_min(rho_min / state.user_to_internal),
            IMobility<dims>(state) {
    }

    ~RegularisedMobility() {}

    void update_mobility() override {
        for(unsigned int idx = 0; idx < this->_sim_state.rho.bins(); idx++) {
            for(int species = 0; species < this->_sim_state.rho.species(); species++) {
                double M = _M * this->_sim_state.rho(idx, species) / (this->_sim_state.rho(idx, species) + _rho_min);
                this->_sim_state.mobility(idx, species) = M;
            }
        }
    }

    GET_NAME(RegularisedMobility)

private:
    double _M;
    double _rho_min;
};

} /* namespace ch */

#endif /* REGULARISEDMOBILITY_H */
