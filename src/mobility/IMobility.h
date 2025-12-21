/*
 * IMobility.h
 *
 * Created on: 12/20/2025
 *      Author: Lorenzo
*/

#ifndef IMOBILITY_H
#define IMOBILITY_H

#include "../SimulationState.h"

namespace ch {

class IMobility {
public:
    IMobility(SimulationState& state) : _sim_state(state) {};

    virtual void update_mobility() = 0;

    virtual ~IMobility() {}

private:
    SimulationState &_sim_state;
};

IMobility *build_mobility(toml::table &config, SimulationState &state);

} /* namespace ch */

#endif /* IMOBILITY_H */
