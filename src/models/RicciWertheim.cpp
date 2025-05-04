/*
 * RicciWertheim.cpp
 *
 *  Created on: May 4, 2025
 *      Author: lorenzo
 */

 #include "RicciWertheim.h"
 #include "../utils/Delta.h"
 
 #ifndef NOCUDA
 #include "../CUDA/models/RicciWertheim.cuh"
 #endif
 
 #include <numeric>
 
 namespace ch {
 
 RicciWertheim::RicciWertheim(toml::table &config) :
				 FreeEnergyModel(config) {
 
	 _B2 = _config_value<double>(config, "ricci.B2");
	 _delta_00 = Delta(config, "ricci.delta_00");
	 _delta_12 = Delta(config, "ricci.delta_12");
 
	 info("B2 = {}, delta_00 = {}, delta_12 = {}", _B2, _delta_00, _delta_12);
 
	 _B2 *= CUB(_user_to_internal);
	 _delta_00 *= CUB(_user_to_internal);
	 _delta_12 *= CUB(_user_to_internal);
 
 #ifndef NOCUDA
	 if(_config_optional_value<bool>(config, "use_CUDA", false)) {
		 init_ricci_symbols(_delta_00, _delta_12);
	 }
 #endif
 }
 
 RicciWertheim::~RicciWertheim() {
 }
 
 double RicciWertheim::bonding_free_energy(const std::vector<double> &rhos) {
	 double X_0 = (-1 + std::sqrt((1 + 12 * rhos[0] * _delta_00))) / (6 * rhos[0] * _delta_00);
	 double A = rhos[1] * _delta_12;
	 double B = rhos[0] * _delta_12;
	 double X_1 = (B - 1 - 2 * A + std::sqrt(SQR(1 + 2 * A - B) + 4 * B)) / (2 * B);
	 double X_2 = (-B + 2 * A - 1 + std::sqrt(SQR(B - 2 * A + 1) + 8 * A)) / (4 * A);

	 double bonding_fe = rhos[0] * (3 * (std::log(X_0) - X_0 / 2) + std::log(X_1) - X_1 / 2 + 2); // nanostars
	 bonding_fe += rhos[1] * (2 * std::log(X_2) - X_2 + 1); // antibodies
 
	 return bonding_fe;
 }
 
 double RicciWertheim::bulk_free_energy(const std::vector<double> &rhos) {
	 double rho = std::accumulate(rhos.begin(), rhos.end(), 0.);
 
	 double mixing_S = 0.;
	 for(int i = 0; i < N_species(); i++) {
		 double x_i = rhos[i] / rho;
		 mixing_S += (x_i > 0) ? x_i * std::log(x_i) : 0;
	 }
	 double B2_contrib = _B2 * rho;
 
	 double f_ref = rho * (std::log(rho * _density_conversion_factor) - 1.0 + mixing_S + B2_contrib);
 
	 return f_ref + bonding_free_energy(rhos);
 }
 
 void RicciWertheim::der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) {
	 for(unsigned int idx = 0; idx < rho.bins(); idx++) {
		 for(int species = 0; species < N_species(); species++) {
			auto rhos = rho.rho_species(idx);
			if(rhos[species] != 0.) {
				double rho_tot = std::accumulate(rhos.begin(), rhos.end(), 0.);
				double der_f_ref = std::log(rhos[species]) + 2.0 * _B2 * rho_tot;

				double der_f_bond;
				double A = rhos[1] * _delta_12;
				double B = rhos[0] * _delta_12;
				if(species == 0) {
					double X_0 = (-1 + std::sqrt((1 + 12 * rhos[0] * _delta_00))) / (6 * rhos[0] * _delta_00);
					double X_1 = (B - 1 - 2 * A + std::sqrt(SQR(1 + 2 * A - B) + 4 * B)) / (2 * B);
					der_f_bond = 3 * std::log(X_0) + std::log(X_1);
				}
				else if(species == 1) {
					double X_2 = (-B + 2 * A - 1 + std::sqrt(SQR(B - 2 * A + 1) + 8 * A)) / (4 * A);
					der_f_bond = 2 * std::log(X_2);
				}
				rho_der(idx, species) = der_f_ref + der_f_bond;
			}
		 }
	 }
 }
 
 void RicciWertheim::der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) {
 #ifndef NOCUDA
	 ricci_wertheim_der_bulk_free_energy(rho, rho_der, vec_size, _B2);
 #endif
 }
 
 } /* namespace ch */
