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

#define SQR(X) ((X) * (X))

struct FreeEnergyModel {
	FreeEnergyModel(cxxopts::ParseResult &options) {

	}

	virtual ~FreeEnergyModel() {

	}

	void print_free_energy(double op_min, double op_max) {
		std::ofstream out("free_energy.dat");

		double step = (op_max - op_min) / 10000;
		for(double op = op_min; op < op_max; op += step) {
			out << op << " " << bulk_free_energy(op) << " " << der_bulk_free_energy(op) << std::endl;
		}

		out.close();
	}

	virtual double der_bulk_free_energy(double op) = 0;
	virtual double bulk_free_energy(double op) = 0;
};

struct LandauModel: public FreeEnergyModel {
	double epsilon;

	LandauModel(cxxopts::ParseResult &options) :
					FreeEnergyModel(options) {
		epsilon = options["epsilon"].as<double>();

		print_free_energy(-1, 1);
	}

	virtual ~LandauModel() {

	}

	double der_bulk_free_energy(double op) override {
		return -epsilon * op + op * op * op;
	}

	double bulk_free_energy(double op) override {
		return -0.5 * epsilon * SQR(op) + 0.25 * SQR(SQR(op));
	}
};

struct WertheimModel: public FreeEnergyModel {
	double B_2 = 2100;
	double valence = 4;
	double T;
	double delta;
	double two_M_delta;

	const double k_B = 1.9872036;

	WertheimModel(cxxopts::ParseResult &options) :
		FreeEnergyModel(options) {
		T = options["temperature"].as<double>();

		double salt = 1.0;
		int L_DNA = 6;

		double delta_H = -42790;
		double delta_S_nosalt = -119.84;
		double delta_S_salt = 0.368 * (L_DNA - 1.0) * std::log(salt);
		double delta_G = delta_H - T * (delta_S_nosalt + delta_S_salt);

		delta = 1.6606 * std::exp(-delta_G / (k_B * T));
		two_M_delta = 2 * valence * delta;

		double psi_average = options["average-psi"].as<double>();
		print_free_energy(psi_average / 10., psi_average * 100);
	}

	virtual ~WertheimModel() {

	}

	double _X(double rho) {
		double sqrt_argument = 2.0 * two_M_delta * rho;
		return (sqrt_argument < 1e-3) ? 1.0 : (-1 + std::sqrt(1 + 2 * two_M_delta * rho)) / (two_M_delta * rho);
	}

	double bulk_free_energy(double rho) override {
		double f_ref = rho * std::log(rho) - rho + B_2 * SQR(rho);
		double f_bond = valence * rho * (std::log(_X(rho)) + 0.5 * (1. - _X(rho)));

		return (f_ref + f_bond);
	}

	double der_bulk_free_energy(double rho) override {
		double der_f_ref = std::log(rho) + 2 * B_2 * rho;
		double X = _X(rho);
		double der_f_bond = valence * (std::log(X) - 0.5 + 1.0 / (2.0 - X) - 0.5 * X / (2.0 - X));

		return (der_f_ref + der_f_bond);
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

	std::vector<double> psi;
	std::vector<std::complex<double>> psi_hat;
	std::vector<double> sqr_wave_vectors, dealiaser;

	fftw_plan psi_plan, psi_inverse_plan;

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

		psi.resize(size);
		psi_hat.resize(size / 2 + 1);
		sqr_wave_vectors.resize(size / 2 + 1);
		dealiaser.resize(size / 2 + 1);

		double psi_average = options["average-psi"].as<double>();
		double psi_noise = options["psi-noise"].as<double>();
		std::generate(psi.begin(), psi.end(), [psi_average, psi_noise]() {
			double noise = 2. * (drand48() - 0.5) * psi_noise;
			return psi_average + noise;
		});

		psi_plan = fftw_plan_dft_r2c_1d(size, psi.data(), reinterpret_cast<fftw_complex *>(psi_hat.data()), FFTW_ESTIMATE);
		// c2r transforms overwrite the input if FFTW_PRESERVE_INPUT is not specified
		psi_inverse_plan = fftw_plan_dft_c2r_1d(size, reinterpret_cast<fftw_complex *>(psi_hat.data()), psi.data(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

		double nyquist_mode = size * M_PI / (N * H) * 2.0 / 3.0;
		for(unsigned int i = 0; i < sqr_wave_vectors.size(); i++) {
			double k = 2.0 * M_PI * i / (N * H);
			dealiaser[i] = (k < nyquist_mode) ? 1.0 : 0.0;
			sqr_wave_vectors[i] = SQR(k);
		}

		fftw_execute(psi_plan);
	}

	virtual ~CahnHilliard() {
		fftw_destroy_plan(psi_plan);
		fftw_destroy_plan(psi_inverse_plan);
	}

	void evolve();
	double total_mass();

	void print_state(std::ofstream &output);
};

template<int dims>
void CahnHilliard<dims>::evolve() {
	static std::vector<double> f_der(size);
	static std::vector<std::complex<double>> f_der_hat(size / 2 + 1);
	static fftw_plan f_plan = fftw_plan_dft_r2c_1d(size, f_der.data(), reinterpret_cast<fftw_complex *>(f_der_hat.data()), FFTW_ESTIMATE);

	for(unsigned int idx = 0; idx < psi.size(); idx++) {
		f_der[idx] = model->der_bulk_free_energy(psi[idx]);
	}

	fftw_execute(f_plan); // transform f_der into f_der_hat

	for(unsigned int k_idx = 0; k_idx < psi_hat.size(); k_idx++) {
		f_der_hat[k_idx] *= dealiaser[k_idx];
		psi_hat[k_idx] = (psi_hat[k_idx] - dt * M * sqr_wave_vectors[k_idx] * f_der_hat[k_idx]) / (1.0 + dt * M * 2.0 * k_laplacian * SQR(sqr_wave_vectors[k_idx]));
	}

	fftw_execute(psi_inverse_plan);

	for(auto &v : psi) {
		v /= size;
	}
}

template<int dims>
double CahnHilliard<dims>::total_mass() {
	return std::accumulate(psi.begin(), psi.end(), 0.);
}

template<>
void CahnHilliard<1>::print_state(std::ofstream &output) {
	for(int idx = 0; idx < size; idx++) {
		output << psi[idx] << " " << std::endl;
	}
	output << std::endl;
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
				modulo <<= bits;
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
	("s,steps", "Number of iterations", cxxopts::value<long long int>())
	("N", "The size of the square grid", cxxopts::value<int>()->default_value("64"))
	("f,free-energy", "The bulk free energy expression to be used (supported values are 'landau' and 'wertheim')", cxxopts::value<std::string>()->default_value("landau"))
	("e,epsilon", "The distance from the critical point in the 'landau' free energy", cxxopts::value<double>()->default_value("0.9"))
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

	FreeEnergyModel *model;
	if(result["free-energy"].as<std::string>() == "landau") {
		model = new LandauModel(result);
	}
	else if(result["free-energy"].as<std::string>() == "wertheim") {
		model = new WertheimModel(result);
	}
	else {
		fprintf(stderr, "ERROR: Unsupported free energy\n");
		exit(1);
	}

	CahnHilliard<1> system(model, result);

	if(result["steps"].count() == 0) {
		fprintf(stderr, "ERROR: The -s/--steps argument in mandatory\n");
		exit(1);
	}
	long long int steps = result["steps"].as<long long int>();
	long long int print_every = result["print-every"].as<long long int>();

	std::ofstream output("init.dat");
	system.print_state(output);
	output.close();

	std::ofstream trajectory;
	if(print_every > 0) {
		trajectory.open("trajectory.dat");
	}
	for(long long int t = 0; t < steps; t++) {
		if(print_every > 0 && t % print_every == 0) {
			system.print_state(trajectory);

			std::ofstream output("last.dat");
			system.print_state(output);
			output.close();

			fprintf(stderr, "%lld %lf %lf\n", t, t * system.dt, system.total_mass());
		}
		system.evolve();
	}
	trajectory.close();

	output.open("last.dat");
	system.print_state(output);
	output.close();

	return 0;
}
