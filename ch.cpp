#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include <cxxopts/cxxopts.hpp>

template<int dims>
struct CahnHilliard {
	int N;
	int N_minus_one;
	int bits;
	int size;
	double dt = 0.01;
	double epsilon = 0.9;
	double k_laplacian = 1.0;
	double D = 1.0;
	double H = 1.0;

	std::vector<double> psi;

	CahnHilliard(int mN, double psi_average) :
					N(mN) {
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

		psi.resize(size);
		std::generate(psi.begin(), psi.end(), [psi_average]() {
			return (drand48() - 0.5) + psi_average;
		});
	}

	void fill_coords(int coords[dims], int idx);
	int cell_idx(int coords[dims]);

	double cell_laplacian(std::vector<double> &field, int idx);
	double bulk_free_energy(double op);

	void evolve();

	void print_state(std::ofstream &output);
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
double CahnHilliard<2>::cell_laplacian(std::vector<double> &field, int idx) {
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

	return (field[cell_idx(coords_xmy)] + field[cell_idx(coords_xpy)] + field[cell_idx(coords_xym)] + field[cell_idx(coords_xyp)] - 4 * field[idx]) / (H * H);
}

template<int dims>
double CahnHilliard<dims>::bulk_free_energy(double op) {
	return epsilon * op - op * op * op;
}

template<int dims>
void CahnHilliard<dims>::evolve() {
	static std::vector<double> free_energy_der(psi.size());
	for(unsigned int idx = 0; idx < psi.size(); idx++) {
		free_energy_der[idx] = k_laplacian * cell_laplacian(psi, idx) + bulk_free_energy(psi[idx]);
	}

	for(unsigned int idx = 0; idx < psi.size(); idx++) {
		psi[idx] += -D * cell_laplacian(free_energy_der, idx) * dt;
	}
}

template<int dims>
void CahnHilliard<dims>::print_state(std::ofstream &output) {
	for(int idx = 0; idx < size; idx++) {
		if(idx > 0) {
			int modulo = N;
			for(int d = 1; d < dims; d++) {
				if(idx % modulo == 0) {
					output << std::endl;
				}
				modulo <<= N;
			}
		}
		output << psi[idx] << " ";
	}
	output << std::endl;
}

int main(int argc, char *argv[]) {
	srand48(51328);

	cxxopts::Options options("cahn-hilliard", "A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation");
	options.add_options()
			("N", "The size of the square grid", cxxopts::value<int>()->default_value("64"))
			("e,epsilon", "The distance from the critical point", cxxopts::value<double>()->default_value("0.9"))
			("dt", "The integration time step", cxxopts::value<double>()->default_value("0.01"))
			("D", "The transport coefficient D of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1.0"))
			("H", "The size of the mesh cells", cxxopts::value<double>()->default_value("1.0"))
			("a,average-psi", "Average value of the order parameter", cxxopts::value<double>()->default_value("0"))
			("s,steps", "Number of iterations", cxxopts::value<long long int>())
			("p,print-every", "Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never)", cxxopts::value<long long int>()->default_value("0"))
			("k", "Strength of the interfacial term of the Cahn-Hilliard equation", cxxopts::value<double>()->default_value("1.0"))
			("h,help", "Print usage")
			;

	auto result = options.parse(argc, argv);

	if(argc == 1 || result.count("help")) {
		fprintf(stderr, "%s", options.help().c_str());
		exit(0);
	}

	CahnHilliard<2> system(result["N"].as<int>(), result["average-psi"].as<double>());
	system.k_laplacian = result["k"].as<double>();
	system.dt = result["dt"].as<double>();
	system.D = result["D"].as<double>();
	system.H = result["H"].as<double>();

	system.epsilon = result["epsilon"].as<double>();

	if(result["steps"].count() == 0) {
		fprintf(stderr, "ERROR: The -s/--steps argument in mandatory\n");
		exit(1);
	}
	long long int steps = result["steps"].as<long long int>();
	long long int print_every = result["print-every"].as<long long int>();

	std::ofstream trajectory;
	if(print_every > 0) {
		trajectory.open("trajectory.dat");
	}
	for(int t = 0; t < steps; t++) {
		if(print_every > 0 && t % print_every == 0) {
			system.print_state(trajectory);
		}
		system.evolve();
	}
	trajectory.close();

	std::ofstream output("last.dat");
	system.print_state(output);
	output.close();

	return 0;
}
