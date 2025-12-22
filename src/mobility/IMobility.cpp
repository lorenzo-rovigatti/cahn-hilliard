/*
 * IMobility.cpp
 *
 * Created on: 20/12/2025
 *     Author: Lorenzo
*/

#include "IMobility.h"

#include "ConstantMobility.h"
namespace ch {

IMobility *build_mobility(toml::table &config, SimulationState &state) {
    std::string mobility_type = config["mobility"]["type"].value_or("constant");

    if(mobility_type == "constant") {
        double M = config["mobility"]["M"].value_or(1.0);
        return new ConstantMobility(state, M);
    }

    throw std::runtime_error("Unknown mobility type: " + mobility_type);
}

} /* namespace ch */
