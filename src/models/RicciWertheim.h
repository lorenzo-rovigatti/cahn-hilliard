/*
 * RicciWertheim.h
 *
 *  Created on: May 4, 2025
 *      Author: lorenzo
 */

 #ifndef SRC_MODELS_RICCIWERTHEIM_H_
 #define SRC_MODELS_RICCIWERTHEIM_H_
 
 #include "FreeEnergyModel.h"
 
 namespace ch {
 
 class RicciWertheim final: public FreeEnergyModel {
 public:
	 RicciWertheim(toml::table &config);
	 virtual ~RicciWertheim();
	 RicciWertheim(const RicciWertheim &other) = default;
	 RicciWertheim(RicciWertheim &&other) = default;
 
	 int N_species() override {
		 return 2;
	 }
 
	 double bonding_free_energy(const std::vector<double> &);
	 double bulk_free_energy(const std::vector<double> &) override;
 
	 void der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) override;
	 void der_bulk_free_energy(const RhoMatrix<double> &rho, RhoMatrix<double> &rho_der) override;
 
	 GET_NAME("Ricci's system Wertheim free energy")
 
 private:
	 double _B2 = 0;
	 double _delta_00, _delta_12;
 };
 
 } /* namespace ch */
 
 #endif /* SRC_MODELS_RICCIWERTHEIM_H_ */
 