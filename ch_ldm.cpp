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

// taken from https://gist.github.com/LingDong-/7e4c4cae5cbbc44400a05fba65f06f23
float ln(float x) {
  unsigned int bx = * (unsigned int *) (&x);
  unsigned int ex = bx >> 23;
  signed int t = (signed int)ex-(signed int)127;
  unsigned int s = (t < 0) ? (-t) : t;
  bx = 1065353216 | (bx & 8388607);
  x = * (float *) (&bx);
  return -1.7417939+(2.8212026+(-1.4699568+(0.44717955-0.056570851*x)*x)*x)*x + 0.6931471806*t;
}

// taken from https://github.com/romeric/fastapprox/blob/master/fastapprox/src/fastlog.h
static inline float fastlog2(float x) {
	union {
		float f;
		uint32_t i;
	} vx = { x };

	union {
		uint32_t i;
		float f;
	} mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };

	float y = vx.i;
	y *= 1.1920928955078125e-7f;

	return y - 124.22551499f
			- 1.498030302f * mx.f
			- 1.72587999f / (0.3520887068f + mx.f);
}

static inline float fastlog(float x) {
	return 0.69314718f * fastlog2(x);
}

#define LR_LOG(x) (std::log((x)))

struct FreeEnergyModel {
	double B_2 = 2100;
	double T;
	double delta_AA, delta_BB;

	const double k_B = 1.9872036;

	FreeEnergyModel(cxxopts::ParseResult &options) {
		T = options["temperature"].as<double>();

		double salt = 1.0;
		int L_DNA = 6;

		double delta_H = -42790;
		double delta_S_nosalt = -119.84;
		double delta_S_salt = 0.368 * (L_DNA - 1.0) * std::log(salt);
		double delta_G = delta_H - T * (delta_S_nosalt + delta_S_salt);

		delta_AA = 1.6606 * std::exp(-delta_G / (k_B * T));
		delta_BB = delta_AA;
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

	double full_numerical_derivative(int species, std::array<double, SPECIES> &partial_rho) {
		double rho_i = partial_rho[species];
		double delta_rho_i = rho_i * 1e-5;

		double fe_r = bulk_free_energy(partial_rho);
		partial_rho[species] += delta_rho_i;
		double fe_rdr = bulk_free_energy(partial_rho);
		partial_rho[species] = rho_i;

		return (fe_rdr - fe_r) / delta_rho_i;
	}

	double der_bulk_free_energy(int species, std::array<double, SPECIES> &partial_rho) {
		// the ideal + B2 part is computed analytically
		double der_f_ref = LR_LOG(partial_rho[species]) + B_2 * partial_rho[species];
		for(int i = 0; i < SPECIES; i++) {
			der_f_ref += B_2 * partial_rho[i];
		}

		// the bonding part is computed numerically
		double rho_i = partial_rho[species];
		double delta_rho_i = rho_i * 1e-5;

		double fe_r = bonding_free_energy(partial_rho);
		partial_rho[species] += delta_rho_i;
		double fe_rdr = bonding_free_energy(partial_rho);
		partial_rho[species] = rho_i;

		return der_f_ref + (fe_rdr - fe_r) / delta_rho_i;
	}

	double bonding_free_energy(std::array<double, SPECIES> &partial_rho) {
		double rho_factor = 1.0 + (2 * partial_rho[2] + partial_rho[4] - 4 * partial_rho[0]) * delta_AA;
		double X_1A = (-rho_factor + std::sqrt(SQR(rho_factor) + 16.0 * partial_rho[0] * delta_AA)) / (8.0 * partial_rho[0] * delta_AA);

		rho_factor = 1.0 + (2 * partial_rho[3] + partial_rho[4] - 4 * partial_rho[1]) * delta_BB;
		double X_2B = (-rho_factor + std::sqrt(SQR(rho_factor) + 16.0 * partial_rho[1] * delta_BB)) / (8.0 * partial_rho[1] * delta_BB);

		double X_3A = 1.0 / (1.0 + 4.0 * partial_rho[0] * X_1A * delta_AA);
		double X_4B = 1.0 / (1.0 + 4.0 * partial_rho[1] * X_2B * delta_BB);

		double f_bond_1 = 4.0 * (LR_LOG(X_1A) + 0.5 * (1. - X_1A));
		double f_bond_2 = 4.0 * (LR_LOG(X_2B) + 0.5 * (1. - X_2B));
		double f_bond_3 = 2.0 * (LR_LOG(X_3A) + 0.5 * (1. - X_3A));
		double f_bond_4 = 2.0 * (LR_LOG(X_4B) + 0.5 * (1. - X_4B));
		double f_bond_5 = 0.5 * (f_bond_3 + f_bond_4);

//		double X_5A = X_3A;
//		double X_5B = X_4B;
//		double f_bond_5 = (std::log(X_5A) - 0.5 * X_5A + std::log(X_5B) - 0.5 * X_5B + 1.0);

		double f_bond =
				partial_rho[0] * f_bond_1 +
				partial_rho[1] * f_bond_2 +
				partial_rho[2] * f_bond_3 +
				partial_rho[3] * f_bond_4 +
				partial_rho[4] * f_bond_5;

		return f_bond;
	}

	double bulk_free_energy(std::array<double, SPECIES> &partial_rho) {
		double rho = std::accumulate(partial_rho.begin(), partial_rho.end(), 0.);

		double mixing_S = 0., B2_contrib = 0.;
		for(int i = 0; i < SPECIES; i++) {
			double x_i = partial_rho[i] / rho;
			mixing_S += x_i * LR_LOG(x_i);

			for(int j = i; j < SPECIES; j++) {
				B2_contrib += B_2 * x_i * partial_rho[j];
			}
		}

		double f_ref = rho * (LR_LOG(rho) - 1.0 + mixing_S + B2_contrib);

		return f_ref + bonding_free_energy(partial_rho);
	}
};

template<int dims>
struct CahnHilliard {
	int N;
	int N_minus_one;
	int bits;
	int size;
	double dt;
	double k_laplacian;
	double M;
	double H;
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

		if(options["load-from"].count() != 0) {
			std::ifstream load_from(options["load-from"].as<std::string>().c_str());

			for(int s = 0; s < SPECIES; s++) {
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
					fprintf(stderr, "Unsupported number of dimensions %d\n", dims);
					exit(1);
				}
			}

			load_from.close();
		}
		else {
			double tetramer_rho = options["tetramer-density"].as<double>();
			double noise = options["noise"].as<double>();
			double R = options["R"].as<double>();
			double linker_rho = 2 * tetramer_rho / (1.0 + R);
			std::for_each(rho.begin(), rho.end(), [this, tetramer_rho, linker_rho, noise, R](std::array<double, SPECIES> &species_rho) {
				species_rho[0] = tetramer_rho * (1 + 2. * (drand48() - 0.5) * noise);
				species_rho[1] = tetramer_rho * (1 + 2. * (drand48() - 0.5) * noise);

				species_rho[2] = linker_rho * (1 + 2. * (drand48() - 0.5) * noise);
				species_rho[3] = linker_rho * (1 + 2. * (drand48() - 0.5) * noise);

				species_rho[4] = 2 * R * linker_rho * (1 + 2. * (drand48() - 0.5) * noise);
			});
		}
	}

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	double cell_laplacian(std::vector<std::array<double, SPECIES>> &field, int species, int idx);

	void evolve();
	double total_mass();

	void print_state(int species, std::ofstream &output);
	void print_density(std::string filename);
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
double CahnHilliard<1>::cell_laplacian(std::vector<std::array<double, SPECIES>> &field, int species, int idx) {
	int idx_m = (idx - 1 + N) & N_minus_one;
	int idx_p = (idx + 1) & N_minus_one;

	return (field[idx_m][species] + field[idx_p][species] - 2.0 * field[idx][species]) / SQR(H);
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
	static std::vector<std::array<double, SPECIES>> rho_der(rho.size());
	// we first evaluate the time derivative for all the fields
	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < SPECIES; species++) {
			rho_der[idx][species] = model->der_bulk_free_energy(species, rho[idx]) - 2 * k_laplacian * cell_laplacian(rho, species, idx);
		}
	}

	// and then we integrate them
	for(unsigned int idx = 0; idx < rho.size(); idx++) {
		for(int species = 0; species < SPECIES; species++) {
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
		output << std::accumulate(rho[idx].begin(), rho[idx].end(), 0.) << " ";
	}
	output << std::endl;

	output.close();
}

int main(int argc, char *argv[]) {
	srand48(51328);

	cxxopts::Options options("cahn-hilliard", "A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation");
	options.add_options()
	("s,steps", "Number of iterations", cxxopts::value<long long int>())
	("N", "The size of the square grid", cxxopts::value<int>()->default_value("64"))
	("l,load-from", "Load the initial conditions from this file", cxxopts::value<std::string>())
	("T,temperature", "Temperature (in Kelvin)", cxxopts::value<double>()->default_value("300"))
	("dt", "The integration time step", cxxopts::value<double>()->default_value("0.001"))
	("M", "The transport coefficient M of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1.0"))
	("H", "The size of the mesh cells", cxxopts::value<double>()->default_value("20.0"))
	("t,tetramer-density", "Average value of the initial density of tetramers", cxxopts::value<double>()->default_value("5e-5"))
	("noise", "Random noise for the initial psi", cxxopts::value<double>()->default_value("0.01"))
	("R", "Fraction of inter-species linkers", cxxopts::value<double>()->default_value("0.1"))
	("p,print-every", "Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never)", cxxopts::value<long long int>()->default_value("0"))
	("k", "Strength of the interfacial term of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1e7"))
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

	for(int i = 0; i < SPECIES; i++) {
		char filename[256];
		sprintf(filename, "init_%d.dat", i);
		std::ofstream output(filename);
		system.print_state(i, output);
		output.close();
	}
	system.print_density("initial_density.dat");

	for(long long int t = 0; t < steps; t++) {
		if(print_every > 0 && t % print_every == 0) {
			for(int i = 0; i < SPECIES; i++) {
				system.print_state(i, trajectory[i]);

				char filename[256];
				sprintf(filename, "last_%d.dat", i);
				std::ofstream output(filename);
				system.print_state(i, output);
				output.close();
			}

			fprintf(stdout, "%lld %lf %lf\n", t, t * system.dt, system.total_mass());
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
	system.print_density("final_density.dat");

	return 0;
}
