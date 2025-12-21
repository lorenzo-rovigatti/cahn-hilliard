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
        state.mobility.fill(M);
    }

    ~ConstantMobility() {}

    void update_mobility() override {
        // Constant mobility: do nothing, since we set it once at the beginning
    }
};

} /* namespace ch */

#endif /* CONSTANTMOBILITY_H */
