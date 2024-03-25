#ifndef SRC_INTEGRATORS_INTEGRATOR_H_
#define SRC_INTEGRATORS_INTEGRATOR_H_

#include "../Object.h"
#include "../models/FreeEnergyModel.h"

namespace ch {

template<int dims>
class Integrator : public Object {
public:
    Integrator(FreeEnergyModel *model, toml::table &config);

    virtual ~Integrator();

    void set_initial_rho(RhoMatrix<double> &rho) override;

    virtual void evolve() = 0;

    GET_NAME(Integrator)

protected:
    int _N_per_dim = 0;
    int _N_bins = 0;
    double _dt = 0.0;
	double _k_laplacian = 0.0;
	double _M = 0.0;
	double _dx = 0.0;
    double _user_to_internal, _internal_to_user;

    FreeEnergyModel *_model;
    RhoMatrix<double> _rho;
};

} /* namespace ch */

#endif