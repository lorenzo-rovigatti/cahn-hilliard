#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>

#include <cxxopts/cxxopts.hpp>

#define SPECIES 5

#define SQR(X) ((X) * (X))

struct FreeEnergyModel {
	double B_2 = 2100;
	double T;
	double delta;

	const double k_B = 1.9872036;

	FreeEnergyModel(cxxopts::ParseResult &options) {
		T = options["temperature"].as<double>();

		double salt = 1.0;
		int L_DNA = 6;

		double delta_H = -42790;
		double delta_S_nosalt = -119.84;
		double delta_S_salt = 0.368 * (L_DNA - 1.0) * std::log(salt);
		double delta_G = delta_H - T * (delta_S_nosalt + delta_S_salt);

		delta = 1.6606 * std::exp(-delta_G / (k_B * T));
	}

	virtual ~FreeEnergyModel() {

	}

	void print_free_energy(double op_min, double op_max) {
		std::ofstream out("free_energy.dat");

		double step = (op_max - op_min) / 10000;
		for(double op = op_min; op < op_max; op += step) {
//			out << op << " " << bulk_free_energy(op) << " " << der_bulk_free_energy(op) << std::endl;
		}

		out.close();
	}

	double der_bulk_free_energy(int species, std::array<double, SPECIES> &partial_rho) {
		return 0.;
	}

	double _X(int species, std::array<double, SPECIES> &partial_rho) {
		return 0.;
	}

	double bulk_free_energy(int species, std::array<double, SPECIES> &partial_rho) {
		double rho = std::accumulate(partial_rho.begin(), partial_rho.end(), 0.);
		double mixing_S = 0., B2_contrib = 0.;
		for(unsigned int i = 0; i < partial_rho.size(); i++) {
			double x_i = partial_rho[i] / rho;
			mixing_S += x_i * std::log(x_i);

			for(unsigned int j = i; j < partial_rho.size(); j++) {
				double x_j = partial_rho[j] / rho;
				B2_contrib += x_i * x_j;
			}
		}
		B2_contrib *= SQR(rho) * B_2;

		double f_ref = rho * std::log(rho) - rho + rho * mixing_S + B2_contrib;

		double rho_factor = 1.0 + (2 * partial_rho[2] + partial_rho[4] - 4 * partial_rho[0]) * delta;
		double X_1A = (-rho_factor + std::sqrt(SQR(rho_factor) + 16.0 * partial_rho[0] * delta)) / (8.0 * partial_rho[0] * delta);

		rho_factor = 1.0 + (2 * partial_rho[3] + partial_rho[5] - 4 * partial_rho[1]) * delta;
		double X_2B = (-(1.0 + rho_factor) + std::sqrt(SQR(1.0 + rho_factor) + 16.0 * partial_rho[1] * delta)) / (8.0 * partial_rho[1] * delta);

		double X_3A = 1.0 / (1.0 + 4.0 * partial_rho[0] * delta);
		double X_4B = 1.0 / (1.0 + 4.0 * partial_rho[1] * delta);

		double X_5A = X_3A;
		double X_5B = X_4B;

		double f_bond =
				partial_rho[0] * 4 * (std::log(X_1A) + 0.5 * (1 - X_1A)) +
				partial_rho[1] * 4 * (std::log(X_2B) + 0.5 * (1 - X_2B)) +
				partial_rho[2] * 2 * (std::log(X_3A) + 0.5 * (1 - X_3A)) +
				partial_rho[3] * 2 * (std::log(X_4B) + 0.5 * (1 - X_4B)) +
				partial_rho[4] * (std::log(X_5A) - 0.5 * X_5A + std::log(X_5B) - 0.5 * X_5B + 1.0);

		return f_ref + f_bond;
	}
};

template<int dims>
struct CahnHilliard {
	int N;
	int N_minus_one;
	int bits;
	int size;
	double dt = 0.01;
	double k_laplacian = 1.0;
	double M = 1.0;
	double H = 1.0;
	FreeEnergyModel *model;

	std::vector<std::array<double, SPECIES>> rho;

	CahnHilliard(FreeEnergyModel *m, cxxopts::ParseResult &options) :
					model(m) {
		N = options["N"].as<int>();
		k_laplacian = options["k"].as<double>();
		dt = options["dt"].as<double>();
		M = options["M"].as<double>();
		H = options["H"].as<double>();

		double log2N = std::log2(N);
		if(ceil(log2N) != floor(log2N)) {
			fprintf(stderr, "N should be a power of 2\n");
			exit(1);
		}

		N_minus_one = N - 1;
		bits = (int) log2N;

		size = N;
		for(int i = 1; i < dims; i++) {
			size *= N;
		}

		rho.resize(size);

		std::for_each(rho.begin(), rho.end(), [this, options](std::array<double, SPECIES> &species_rho) {
			double psi_average = options["average-psi"].as<double>();
			double psi_noise = options["psi-noise"].as<double>();
			std::generate(species_rho.begin(), species_rho.end(), [psi_average, psi_noise]() {
				double noise = 2. * (drand48() - 0.5) * psi_noise;
				return psi_average + noise;
			});
		});

	}

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	double cell_laplacian(std::vector<std::array<double, SPECIES>> &field, int species, int idx);

	void evolve();
	double total_mass();

	void print_state(int species, std::ofstream &output);
};

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
double CahnHilliard<2>::cell_laplacian(std::vector<std::array<double, SPECIES>> &field, int species, int idx) {
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
	static std::vector<std::array<double, SPECIES>> free_energy_der(rho.size());
	for(int species = 0; species < SPECIES; species++) {
		for(unsigned int idx = 0; idx < rho[species].size(); idx++) {
			free_energy_der[idx][species] = model->der_bulk_free_energy(species, rho[idx]) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
		}

		for(unsigned int idx = 0; idx < rho[species].size(); idx++) {
			rho[idx][species] += M * cell_laplacian(free_energy_der, species, idx) * dt;
		}
	}
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	double mass = 0.;
	for(int i = 0; i < rho.size(); i++) {
		mass += std::accumulate(rho[i].begin(), rho[i].end(), 0.);
	}

	return mass;
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
		output << rho[species][idx] << " ";
	}
	output << std::endl;
}

int main(int argc, char *argv[]) {
	srand48(51328);

	cxxopts::Options options("cahn-hilliard", "A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation");
	options.add_options()
	("s,steps", "Number of iterations", cxxopts::value<long long int>())
	("N", "The size of the square grid", cxxopts::value<int>()->default_value("64"))
	("T,temperature", "Temperature (in Kelvin), used by the 'wertheim' free energy", cxxopts::value<double>()->default_value("300"))
	("dt", "The integration time step", cxxopts::value<double>()->default_value("0.01"))
	("M", "The transport coefficient M of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1.0"))
	("H", "The size of the mesh cells", cxxopts::value<double>()->default_value("1.0"))
	("a,average-psi", "Average value of the order parameter", cxxopts::value<double>()->default_value("0"))
	("psi-noise", "Random noise for the initial psi", cxxopts::value<double>()->default_value("0.1"))
	("p,print-every", "Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never)", cxxopts::value<long long int>()->default_value("0"))
	("k", "Strength of the interfacial term of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1.0"))
	("h,help", "Print usage");

	auto result = options.parse(argc, argv);

	if(argc == 1 || result.count("help")) {
		fprintf(stderr, "%s", options.help().c_str());
		exit(0);
	}

	FreeEnergyModel *model = new FreeEnergyModel(result);
	CahnHilliard<2> system(model, result);

	if(result["steps"].count() == 0) {
		fprintf(stderr, "ERROR: The -s/--steps argument in mandatory\n");
		exit(1);
	}
	long long int steps = result["steps"].as<long long int>();
	long long int print_every = result["print-every"].as<long long int>();

	std::ofstream trajectory[SPECIES];
	if(print_every > 0) {
		for(int i = 0; i < SPECIES; i++) {
			char filename[256];
			sprintf(filename, "trajectory_%d.dat", i);
			trajectory[i].open(filename);
		}
	}

	for(int t = 0; t < steps; t++) {
		if(print_every > 0 && t % print_every == 0) {
			for(int i = 0; i < SPECIES; i++) {
				system.print_state(i, trajectory[i]);
			}

			fprintf(stdout, "%d %lf %lf\n", t, t * system.dt, system.total_mass());
		}
		system.evolve();
	}

	if(print_every > 0) {
		for(int i = 0; i < SPECIES; i++) {
			trajectory[i].close();
		}
	}

	for(int i = 0; i < SPECIES; i++) {
		char filename[256];
		sprintf(filename, "last_%d.dat", i);
		std::ofstream output(filename);
		system.print_state(i, output);
		output.close();
	}

	return 0;
}
