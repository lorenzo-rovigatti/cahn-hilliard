#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>

#include <complex>
#include <fftw3.h>

#include <cxxopts/cxxopts.hpp>

#define SPECIES 2

#define SQR(X) ((X) * (X))
#define LR_LOG(x) (std::log((x)))

struct FreeEnergyModel {
	double B_2 = 2100 * 0.5;
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

		fprintf(stderr, "delta_AA = %lf, delta_BB = %lf\n", delta_AA, delta_BB);
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
		double der_f_ref = LR_LOG(partial_rho[species]);
		for(int i = 0; i < SPECIES; i++) {
			der_f_ref += 2 * B_2 * partial_rho[i];
		}

		// the bonding part is computed numerically
		double rho_i = partial_rho[species];
		double delta_rho_i = rho_i * 1e-4;

		double fe_r = bonding_free_energy(partial_rho);
		partial_rho[species] += delta_rho_i;
		double fe_rdr = bonding_free_energy(partial_rho);
		partial_rho[species] = rho_i;

		return der_f_ref + (fe_rdr - fe_r) / delta_rho_i;
	}

	double bonding_free_energy(std::array<double, SPECIES> &partial_rho) {
		double AA = 4.0 * partial_rho[0] * delta_AA;
		double BB = 1.0 + (2 * partial_rho[1] - 4.0 * partial_rho[0]) * delta_AA;
		double CC = -1.0;
		double X_1A = (-BB + std::sqrt(SQR(BB) - 4.0 * AA * CC)) / (2.0 * AA);

		double X_3A = 1.0 / (1.0 + 4.0 * partial_rho[0] * X_1A * delta_AA);

		double f_bond_1 = 4.0 * (LR_LOG(X_1A) + 0.5 * (1. - X_1A));
		double f_bond_2 = 2.0 * (LR_LOG(X_3A) + 0.5 * (1. - X_3A));

		double f_bond =
				partial_rho[0] * f_bond_1 +
				partial_rho[1] * f_bond_2;

		return f_bond;
	}

	double bulk_free_energy(std::array<double, SPECIES> &partial_rho) {
		double rho = std::accumulate(partial_rho.begin(), partial_rho.end(), 0.);

		double mixing_S = 0., B2_contrib = 0.;
		for(int i = 0; i < SPECIES; i++) {
			double x_i = partial_rho[i] / rho;
			mixing_S += partial_rho[i] * LR_LOG(x_i);

			for(int j = 0; j < SPECIES; j++) {
				B2_contrib += B_2 * partial_rho[i] * partial_rho[j];
			}
		}

		double f_ref = rho * LR_LOG(rho) - rho + mixing_S + B2_contrib;

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

	std::array<std::vector<double>, SPECIES> rho;
	std::array<std::vector<std::complex<double>>, SPECIES> rho_hat;
	std::vector<double> sqr_wave_vectors, dealiaser;

	std::array<fftw_plan, SPECIES> rho_plans, rho_inverse_plans;

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

		for(int species = 0; species < SPECIES; species++) {
			rho[species].resize(size);
		}

		sqr_wave_vectors.resize(size / 2 + 1);
		dealiaser.resize(size / 2 + 1);

		if(options["load-from"].count() != 0) {
			std::ifstream load_from(options["load-from"].as<std::string>().c_str());

			for(int s = 0; s < SPECIES; s++) {
				switch(dims) {
				case 1:
					for(int idx = 0; idx < N; idx++) {
						load_from >> rho[s][idx];
					}
					break;
				case 2:
					int coords[2];
					for(coords[1] = 0; coords[1] < N; coords[1]++) {
						for(coords[0] = 0; coords[0] < N; coords[0]++) {
							int idx = cell_idx(coords);
							load_from >> rho[s][idx];
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
			double linker_rho = 2 * tetramer_rho;

			for(int idx = 0; idx < size; idx++) {
				rho[0][idx] = tetramer_rho * (1 + 2. * (drand48() - 0.5) * noise);
				rho[1][idx] = linker_rho * (1 + 2. * (drand48() - 0.5) * noise);
			}
		}

		for(int species = 0; species < SPECIES; species++) {
			rho_hat[species].resize(size / 2 + 1);
			rho_plans[species] = fftw_plan_dft_r2c_1d(size, rho[species].data(), reinterpret_cast<fftw_complex *>(rho_hat[species].data()), FFTW_ESTIMATE);
			// c2r transforms overwrite the input if FFTW_PRESERVE_INPUT is not specified
			rho_inverse_plans[species] = fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex *>(rho_hat[species].data()), rho[species].data(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

			fftw_execute(rho_plans[species]);
		}

		double nyquist_mode = size * M_PI / (N * H) * 2.0 / 3.0;
		for(unsigned int i = 0; i < sqr_wave_vectors.size(); i++) {
			double k = 2.0 * M_PI * i / (N * H);
			dealiaser[i] = (k < nyquist_mode) ? 1.0 : 0.0;
			sqr_wave_vectors[i] = SQR(k);
		}
	}

	int cell_idx(int coords[dims]);

	void evolve();
	double total_mass();

	void print_state(int species, std::ofstream &output);
	void print_density(std::string filename);
};

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

template<int dims>
void CahnHilliard<dims>::evolve() {
	static std::array<std::vector<double>, SPECIES> f_der = {
			std::vector<double>(size),
			std::vector<double>(size),
	};

	static std::array<std::vector<std::complex<double>>, SPECIES> f_der_hat = {
			std::vector<std::complex<double>>(size / 2 + 1),
			std::vector<std::complex<double>>(size / 2 + 1)
	};

	static std::array<fftw_plan, SPECIES> f_plans = {
			fftw_plan_dft_r2c_1d(size, f_der[0].data(), reinterpret_cast<fftw_complex *>(f_der_hat[0].data()), FFTW_ESTIMATE),
			fftw_plan_dft_r2c_1d(size, f_der[1].data(), reinterpret_cast<fftw_complex *>(f_der_hat[1].data()), FFTW_ESTIMATE)
	};

	for(int species = 0; species < SPECIES; species++) {
		for(unsigned int idx = 0; idx < rho[species].size(); idx++) {
			std::array<double, SPECIES> partial_rho = {
				rho[0][idx],
				rho[1][idx]
			};
			f_der[species][idx] = model->der_bulk_free_energy(species, partial_rho);
		}
		fftw_execute(f_plans[species]);
	}

	for(int species = 0; species < SPECIES; species++) {
		for(unsigned int k_idx = 0; k_idx < rho_hat[species].size(); k_idx++) {
			f_der_hat[species][k_idx] *= dealiaser[k_idx];
			rho_hat[species][k_idx] = (rho_hat[species][k_idx] - dt * M * sqr_wave_vectors[k_idx] * f_der_hat[species][k_idx]) / (1.0 + dt * M * 2.0 * k_laplacian * SQR(sqr_wave_vectors[k_idx]));
		}

		fftw_execute(rho_inverse_plans[species]);

		for(auto &v : rho[species]) {
			v /= size;
		}
	}
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	double mass = 0.;
	for(int idx = 0; idx < size; idx++) {
		mass += rho[0][idx] + rho[1][idx];
	}

	return mass;
}

template<>
void CahnHilliard<1>::print_state(int species, std::ofstream &output) {
	for(int idx = 0; idx < size; idx++) {
		output << rho[species][idx] << " " << std::endl;
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
		output << rho[species][idx] << " ";
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
		output << rho[0][idx] + rho[1][idx] << " ";
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
	("p,print-every", "Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never)", cxxopts::value<long long int>()->default_value("0"))
	("k", "Strength of the interfacial term of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1e7"))
	("h,help", "Print usage");

	auto result = options.parse(argc, argv);

	if(argc == 1 || result.count("help")) {
		fprintf(stderr, "%s", options.help().c_str());
		exit(0);
	}

	std::ofstream cmd_output("executed_command");
	for(int i = 0; i < argc; i++) {
		cmd_output << argv[i] << " ";
	}
	cmd_output << std::endl;
	cmd_output.close();

	FreeEnergyModel *model = new FreeEnergyModel(result);
	CahnHilliard<1> system(model, result);

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
