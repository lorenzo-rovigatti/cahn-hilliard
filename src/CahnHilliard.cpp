/*
 * CahnHilliard.cpp
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#include "CahnHilliard.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace ch {

template<int dims>
CahnHilliard<dims>::CahnHilliard(FreeEnergyModel *m, toml::table &config) :
				model(m) {

	N = config["N"].value_or(64);
	k_laplacian = config["k"].value_or(1.0);
	dt = config["dt"].value_or(0.001);
	M = config["M"].value_or(1.0);
	H = config["H"].value_or(1.0);

	double log2N = std::log2(N);
	if(ceil(log2N) != floor(log2N)) {
		critical("N should be a power of 2");
	}

	N_minus_one = N - 1;
	bits = (int) log2N;

	size = N;
	for(int i = 1; i < dims; i++) {
		size *= N;
	}

	rho.resize(size, std::vector<double>(model->N_species(), 0.));

	if(!config["initial-density"] && !config["load-from"]) {
		critical("Either 'initial-density' or 'load-from' should be specified");
	}

	if(config["load-from"]) {
		std::ifstream load_from(config["load-from"].value<std::string>().value().c_str());

		for(int s = 0; s < model->N_species(); s++) {
			switch(dims) {
			case 1:
				for(int idx = 0; idx < N; idx++) {
					load_from >> rho[idx][s];
				}
				break;
			case 2:
				int coords[2];
				for(coords[1] = 0; coords[1] < N; coords[1]++) {
					for(coords[0] = 0; coords[0] < N; coords[0]++) {
						int idx = cell_idx(coords);
						load_from >> rho[idx][s];
					}
				}

				break;
			default:
				critical("Unsupported number of dimensions {}", dims);
			}
		}

		load_from.close();
	}
	else { // initial-density
		auto densities = config["initial-density"];
		if(model->N_species() > 1) {
			if(!densities.is_array() || densities.as_array()->size()) {
				critical("initial-density should contain as many elements as the number of species ({})", model->N_species());
			}

			toml::array& rho_array = *densities.as_array();

			std::for_each(rho.begin(), rho.end(), [this, rho_array](std::vector<double> &species_rho) {
				for(int i = 0; i < species_rho.size(); i++) {
					double average_rho = *rho_array.get_as<double>(i);
					species_rho[i] = average_rho * (1.0 + 2.0 * (drand48() - 0.5) * 1e-2);
				}
			});
		}
		else { // single species
			double average_rho = *densities.as_floating_point();
			std::for_each(rho.begin(), rho.end(), [this, average_rho](std::vector<double> &species_rho) {
				species_rho[0] = average_rho * (1.0 + 2.0 * (drand48() - 0.5) * 1e-2);
			});
		}
	}
}

template<int dims>
CahnHilliard<dims>::~CahnHilliard() {

}

template<int dims>
void CahnHilliard<dims>::fill_coords(int coords[dims], int idx) {
	for(int d = 0; d < dims; d++) {
		coords[d] = idx & N_minus_one;
		idx >>= bits; // divide by N
	}
}

template<int dims>
int CahnHilliard<dims>::cell_idx(int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= bits; // multiply by N
	}
	return idx;
}

template<>
double CahnHilliard<1>::cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx) {
	int idx_m = (idx - 1 + N) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (field[idx_m][species] + field[idx_p][species] - 2.0 * field[idx][species]) / SQR(H);
}

template<>
double CahnHilliard<2>::cell_laplacian(std::vector<std::vector<double>> &field, int species, int idx) {
	int coords_xy[2];
	fill_coords(coords_xy, idx);

	int coords_xmy[2] = {
			(coords_xy[0] - 1 + N) & N_minus_one,
			coords_xy[1]
	};

	int coords_xym[2] = {
			coords_xy[0],
			(coords_xy[1] - 1 + N) & N_minus_one
	};

	int coords_xpy[2] = {
			(coords_xy[0] + 1) & N_minus_one,
			coords_xy[1]
	};

	int coords_xyp[2] = {
			coords_xy[0],
			(coords_xy[1] + 1) & N_minus_one
	};

	return (
			field[cell_idx(coords_xmy)][species] +
			field[cell_idx(coords_xpy)][species] +
			field[cell_idx(coords_xym)][species] +
			field[cell_idx(coords_xyp)][species] -
			4 * field[idx][species])
			/ SQR(H);
}

template<int dims>
void CahnHilliard<dims>::evolve() {
	static std::vector<std::vector<double>> rho_der(rho.size(), std::vector<double>(model->N_species()));
	// we first evaluate the time derivative for all the fields
	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho_der[idx][species] = model->der_bulk_free_energy(species, rho[idx]) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
		}
	}

	// and then we integrate them
	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < model->N_species(); species++) {
			rho[idx][species] += M * cell_laplacian(rho_der, species, idx) * dt;
		}
	}
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	double mass = 0.;
	for(unsigned int i = 0; i < rho.size(); i++) {
		mass += std::accumulate(rho[i].begin(), rho[i].end(), 0.);
	}

	return mass;
}

template<>
void CahnHilliard<1>::print_state(int species, std::ofstream &output) {
	for(int idx = 0; idx < size; idx++) {
		output << rho[idx][species] << " " << std::endl;
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_state(int species, std::ofstream &output) {
	for(int idx = 0; idx < size; idx++) {
		if(idx > 0) {
			int modulo = N;
			for(int d = 1; d < dims; d++) {
				if(idx % modulo == 0) {
					output << std::endl;
				}
				modulo <<= bits;
			}
		}
		output << rho[idx][species] << " ";
	}
	output << std::endl;
}

template<int dims>
void CahnHilliard<dims>::print_density(std::string filename) {
	std::ofstream output(filename.c_str());

	for(int idx = 0; idx < size; idx++) {
		if(idx > 0) {
			int modulo = N;
			for(int d = 1; d < dims; d++) {
				if(idx % modulo == 0) {
					output << std::endl;
				}
				modulo <<= bits;
			}
		}
		output << std::accumulate(rho[idx].begin(), rho[idx].end(), 0.) << std::endl;
	}

	output.close();
}

template class CahnHilliard<1>;
template class CahnHilliard<2>;

} /* namespace ch */
