/*
 * ConstantMobility.h
 *
 * Created on: 12/20/2025
 *      Author: Lorenzo
*/

#ifndef CONSTANTMOBILITY_H
#define CONSTANTMOBILITY_H

#include "IMobility.h"

namespace ch {

class ConstantMobility: public IMobility {
public:
    ConstantMobility(SimulationState& state, double M) : 
            IMobility(state) {
        M /= state.user_to_internal; // proportional to m^-1
        state.mobility.fill(M);
    }

    ~ConstantMobility() {}

    void update_mobility() override {
        // we don't need to do anything, since the mobility is set once at the beginning
    }

    GET_NAME(ConstantMobility)
};

} /* namespace ch */

#endif /* CONSTANTMOBILITY_H */
