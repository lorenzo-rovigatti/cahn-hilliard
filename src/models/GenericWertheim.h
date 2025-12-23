/*
 * GenericWertheim.h
 *
 *  Created on: Feb 19, 2025
 *      Author: lorenzo
 * 
 *  This file implements a generic multi-species Wertheim free-energy model
 *  for associating (patchy) particles.
 *
 *  Each species is characterized by a set of bonding patches. Patches are
 *  identified by integer labels, and bonding between patches is controlled
 *  by a matrix of bonding volumes (delta).
 *
 *  The implementation supports:
 *   - multiple species
 *   - arbitrary patch valence per species
 *   - patch–patch bonding rules
 *   - species-dependent virial coefficients (B2)
 */

#ifndef SRC_MODELS_GENERICWERTHEIM_H_
#define SRC_MODELS_GENERICWERTHEIM_H_

#include "FreeEnergyModel.h"

#include <set>

namespace ch {

/**
 * @brief Representation of a unique patch type within a species.
 *
 * A species may contain multiple identical patches (same patch index).
 * This struct stores:
 *  - idx:          the global patch identifier
 *  - multiplicity: how many times this patch appears on the species
 *
 * Example:
 *   patches = [0, 0, 1]
 *   unique_patches = [{idx=0, multiplicity=2}, {idx=1, multiplicity=1}]
 */
struct UniquePatch {
	int idx;
	int multiplicity;
};

/**
 * @brief Patch-level interaction with another species.
 *
 * This structure is used when iterating over the bonding partners of a
 * given patch type.
 *
 * For a given patch A, a PatchInteraction stores:
 *  - species: index of the interacting species
 *  - patches: list of unique patches on that species that can bond to A,
 *             together with their multiplicities
 *
 * This allows efficient evaluation of the Wertheim mass-action equations.
 */
struct PatchInteraction {
	int species = -1;
	std::vector<UniquePatch> patches;
};

/**
 * @brief Description of a particle species.
 *
 * Each species is defined by:
 *  - a unique species index
 *  - a list of all patches carried by the species (including repetitions)
 *  - a compressed list of unique patches with multiplicities
 *
 * The patch indices refer to a global patch labeling shared by all species.
 */
struct Species {
	int idx;
	int N_unique_patches;

	// list of all patches on the species (may contain repetitions),
	// e.g. [0, 0, 1] for a particle with two patches of type 0 and one of type 1
	std::vector<int> patches;

	// list of unique patches with multiplicities,
	// e.g. [{0,2}, {1,1}] for the example above
	std::vector<UniquePatch> unique_patches;
};

/**
 * @brief Term in the derivative of the free energy with respect to density.
 * 
 * This struct is used to precompute contributions to the fraction of unbonded
 * patches when calculating the free energy and its derivative.
 */
struct Term { 
	int species; 
	int other_patch; 
	double coeff; 
}; 


/**
 * @brief Generic multi-species Wertheim free-energy model.
 *
 * This class implements the bulk free energy and its derivatives for a
 * system of associating particles using Wertheim's first-order thermodynamic
 * perturbation theory.
 *
 * The model supports:
 *  - multiple species with different patch compositions
 *  - arbitrary patch–patch bonding rules
 *  - species-dependent second virial coefficients
 *
 * The class also provides a mobility consistent with diffusive dynamics
 * of conserved densities.
 */
class GenericWertheim final: public FreeEnergyModel {
public:
	GenericWertheim(toml::table &config);
	virtual ~GenericWertheim();
	GenericWertheim(const GenericWertheim &other) = default;
	GenericWertheim(GenericWertheim &&other) = default;

	int N_species() override {
		return _species.size();
	}

	double bonding_free_energy(const SpeciesView<double> &);
	double bulk_free_energy(const SpeciesView<double> &) override;

	void der_bulk_free_energy(field_type *rho, float *rho_der, int vec_size) override;
	void der_bulk_free_energy(const MultiField<double> &rho, MultiField<double> &rho_der) override;

	void set_mobility(const MultiField<double> &rho, double M0, MultiField<double> &mobility) override;

	GET_NAME("Generic Wertheim free energy")

private:
	/// List of species in the system
	std::vector<Species> _species;
	/// Global list of unique patch identifiers appearing in any species
	std::vector<int> _unique_patch_ids;
	/**
	 * For each unique patch, a list of PatchInteraction objects describing
	 * which species and which patches can bond to it.
	 *
	 * Indexed as: _unique_patch_interactions[patch_idx]
	 */
	std::vector<std::vector<PatchInteraction>> _unique_patch_interactions;
	/**
	 * Precomputed terms for calculations of the fraction of unbonded patches
	 *
	 * Indexed as: _terms_per_patch[patch_idx] = list of Term
	 */
	std::vector<std::vector<Term>> _terms_per_patch;
	/**
	 * Bonding volume matrix (delta):
	 *   delta[patch_A * N_patches + patch_B]
	 *
	 * Symmetric, zero if patches do not interact.
	 */
	std::vector<double> _delta;
	/**
	 * Second virial coefficient matrix (B2):
	 *   B2[species_A * N_species + species_B]
	 *
	 * Symmetric in species indices.
	 */
	std::vector<double> _B2;
	/// Total number of distinct patch types in the system
	int _N_patches = 0;
	/// Preallocated storage for mass-action equation solutions
	std::vector<std::vector<double>> _all_Xs;

	std::pair<int, int> _parse_interaction(std::string int_string, std::string context);
	void _update_X(const SpeciesView<double> &, std::vector<double> &);
	double _der_contribution(const std::vector<double> &, int);
};

} /* namespace ch */

#endif /* SRC_MODELS_GENERICWERTHEIM_H_ */
