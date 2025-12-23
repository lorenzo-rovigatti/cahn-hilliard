/*
 * IMobility.cpp
 *
 * Created on: 20/12/2025
 *     Author: Lorenzo
*/

#include "IMobility.h"

#include "ConstantMobility.h"
#include "FreeEnergyMobility.h"
#include "RegularisedMobility.h"

namespace ch {

IMobility *build_mobility(toml::table &config, SimulationState &state) {
    const std::string error_source = "build_mobility";
    std::string mobility_type = config["mobility"]["type"].value_or("constant");
    double M = config["mobility"]["M"].value_or(1.0);

    if(mobility_type == "constant") {
        return new ConstantMobility(state, M);
    }
    else if(mobility_type == "regularised") {
        auto rho_opt = config["mobility"]["rho_min"].value<double>();
        if(rho_opt == std::nullopt) {
            spdlog::critical("For 'regularised' mobility, 'rho_min' must be specified (error source: {})", error_source);
            exit(1);
        }
        double rho_min = rho_opt.value();
        return new RegularisedMobility(state, M, rho_min);
    }
    else if(mobility_type == "free_energy") {
        return new FreeEnergyMobility(state, M);
    }

    throw std::runtime_error("Unknown mobility type: " + mobility_type);
}

} /* namespace ch */
