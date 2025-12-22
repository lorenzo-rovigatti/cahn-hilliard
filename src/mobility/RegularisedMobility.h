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

class RegularisedMobility: public IMobility {
public:
    RegularisedMobility(SimulationState& state, double M, double rho_min) :
            _M(M / state.user_to_internal), 
            _rho_min(rho_min / state.user_to_internal),
            IMobility(state) {
    }

    ~RegularisedMobility() {}

    void update_mobility() override {
        for(unsigned int idx = 0; idx < _sim_state.rho.bins(); idx++) {
            for(int species = 0; species < _sim_state.rho.species(); species++) {
                double M = _M * _sim_state.rho(idx, species) / (_sim_state.rho(idx, species) + _rho_min);
                _sim_state.mobility(idx, species) = M;
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
