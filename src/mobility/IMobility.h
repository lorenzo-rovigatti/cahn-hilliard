/*
 * IMobility.h
 *
 * Created on: 12/20/2025
 *      Author: Lorenzo
*/

#ifndef IMOBILITY_H
#define IMOBILITY_H

#include "../defs.h"
#include "../Object.h"
#include "../SimulationState.h"

namespace ch {

template<int dims>
class IMobility : public Object {
public:
    IMobility(SimulationState<dims>& state) : _sim_state(state) {};

    virtual void update_mobility() = 0;

    virtual ~IMobility() {}

    GET_NAME(IMobility)

protected:
    SimulationState<dims> &_sim_state;
};

// Builders for supported dimensions
IMobility<1> *build_mobility(toml::table &config, SimulationState<1> &state);
IMobility<2> *build_mobility(toml::table &config, SimulationState<2> &state);

} /* namespace ch */

#endif /* IMOBILITY_H */
