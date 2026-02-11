/*
 * IMobility.cpp
 *
 * Created on: 20/12/2025
 *     Author: Lorenzo
*/

#include "IMobility.h"

#include "ConstantMobility.h"
#include "FreeEnergyMobility.h"
#include "GelMobility.h"
#include "RegularisedMobility.h"

namespace ch {

IMobility<1> *build_mobility(toml::table &config, SimulationState<1> &state) {
    const std::string error_source = "build_mobility";
    std::string mobility_type = config["mobility"]["type"].value_or("constant");
    double M = config["mobility"]["M"].value_or(1.0);

    if(mobility_type == "constant") {
        return new ConstantMobility<1>(state, M);
    }
    else if(mobility_type == "free_energy") {
        return new FreeEnergyMobility<1>(state, M);
    }
    else if(mobility_type == "gel") {
        double dt = config["dt"].value<double>().value();
        auto phi_critical_opt = config["mobility"]["phi_critical"].value<double>();
        auto c_0_opt = config["mobility"]["c_0"].value<double>();
        auto M_c_opt = config["mobility"]["M_c"].value<double>();
        auto epsilon_opt = config["landau"]["epsilon"].value<double>();
        if(phi_critical_opt == std::nullopt || c_0_opt == std::nullopt || M_c_opt == std::nullopt || epsilon_opt == std::nullopt) {
            spdlog::critical("For 'gel' mobility, 'phi_critical', 'c_0', 'M_c', and 'landau.epsilon' must be specified (error source: {})", error_source);
            exit(1);
        }
        double epsilon = epsilon_opt.value();
        double beta_delta_F = 10.0 / (1 - epsilon) - std::log(24000.0);
        double phi_critical = phi_critical_opt.value();
        double c_0 = c_0_opt.value();
        double M_c = M_c_opt.value();
        return new GelMobility<1>(state, dt, phi_critical, c_0, M_c, beta_delta_F);
    }
    else if(mobility_type == "regularised") {
        auto rho_opt = config["mobility"]["rho_min"].value<double>();
        if(rho_opt == std::nullopt) {
            spdlog::critical("For 'regularised' mobility, 'rho_min' must be specified (error source: {})", error_source);
            exit(1);
        }
        double rho_min = rho_opt.value();
        return new RegularisedMobility<1>(state, M, rho_min);
    }

    throw std::runtime_error("Unknown mobility type: " + mobility_type);
}

IMobility<2> *build_mobility(toml::table &config, SimulationState<2> &state) {
    const std::string error_source = "build_mobility";
    std::string mobility_type = config["mobility"]["type"].value_or("constant");
    double M = config["mobility"]["M"].value_or(1.0);

    if(mobility_type == "constant") {
        return new ConstantMobility<2>(state, M);
    }
    else if(mobility_type == "free_energy") {
        return new FreeEnergyMobility<2>(state, M);
    }
    else if(mobility_type == "gel") {
        double dt = config["dt"].value<double>().value();
        auto phi_critical_opt = config["mobility"]["phi_critical"].value<double>();
        auto c_0_opt = config["mobility"]["c_0"].value<double>();
        auto M_c_opt = config["mobility"]["M_c"].value<double>();
        auto epsilon_opt = config["landau"]["epsilon"].value<double>();
        if(phi_critical_opt == std::nullopt || c_0_opt == std::nullopt || M_c_opt == std::nullopt || epsilon_opt == std::nullopt) {
            spdlog::critical("For 'gel' mobility, 'phi_critical', 'c_0', 'M_c', and 'landau.epsilon' must be specified (error source: {})", error_source);
            exit(1);
        }
        double epsilon = epsilon_opt.value();
        double beta_delta_F = 10.0 / (1 - epsilon) - std::log(24000.0);
        double phi_critical = phi_critical_opt.value();
        double c_0 = c_0_opt.value();
        double M_c = M_c_opt.value();
        return new GelMobility<2>(state, dt, phi_critical, c_0, M_c, beta_delta_F);
    }
    else if(mobility_type == "regularised") {
        auto rho_opt = config["mobility"]["rho_min"].value<double>();
        if(rho_opt == std::nullopt) {
            spdlog::critical("For 'regularised' mobility, 'rho_min' must be specified (error source: {})", error_source);
            exit(1);
        }
        double rho_min = rho_opt.value();
        return new RegularisedMobility<2>(state, M, rho_min);
    }

    throw std::runtime_error("Unknown mobility type: " + mobility_type);
}

} /* namespace ch */
