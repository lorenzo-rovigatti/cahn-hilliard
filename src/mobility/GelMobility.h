/*
 * GelMobility.h
 *
 * Created on: 12/22/2025
 *      Author: Lorenzo
*/

#ifndef GELMOBILITY_H
#define GELMOBILITY_H

#include "IMobility.h"

#include <random>

namespace ch {

template<int dims>
class GelMobility: public IMobility<dims> {
public:
    GelMobility(SimulationState<dims>& state, double dt, double phi_critical, double c_0, double M_c, double beta_delta_F) :
            _dt(dt),
            _phi_critical(phi_critical),
            _c_0(c_0),
            _M_c(M_c / state.user_to_internal),
            IMobility<dims>(state) {

        _p_gel = exp(beta_delta_F) / (1.0 + exp(beta_delta_F));

        _gel_OP = MultiField<double>(state.mobility.bins(), 1);
        std::mt19937 generator;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for(unsigned int idx = 0; idx < state.mobility.bins(); idx++) {
            _gel_OP(idx, 0) = dist(generator) * 1e-4;
        }

        this->info("p_gel = {}, phi_critical = {}, c_0 = {}, M_c = {}, beta_delta_F = {}", _p_gel, _phi_critical, _c_0, _M_c, beta_delta_F);
    }

    ~GelMobility() {}

    void update_mobility() override {
        for(unsigned int idx = 0; idx < this->_sim_state.rho.bins(); idx++) {
            double rho_tot = 0.;
            for(int species = 0; species < this->_sim_state.rho.species(); species++) {
                rho_tot += this->_sim_state.rho(idx, species);
            }
            
            // update the gel OP
            double c = _gel_OP(idx, 0);
            double phi = (rho_tot + 1.0) / 2.0;
            double g = (_p_gel * phi - _phi_critical) / (1.0 - _phi_critical);
            double c_der = _M_c * (g * c - c * c);

            _gel_OP(idx, 0) += c_der * this->_dt;
            for(int species = 0; species < this->_sim_state.rho.species(); species++) {
                double M = std::exp(-_gel_OP(idx, 0) / _c_0);
                this->_sim_state.mobility(idx, species) = M;
            }
        }
    }

    GET_NAME(GelMobility)

private:
    double _dt;
    double _phi_critical;
    double _c_0;
    double _M_c;
    double _p_gel;
    double _beta_delta_F;

    MultiField<double> _gel_OP;
};

} /* namespace ch */

#endif /* GELMOBILITY_H */
