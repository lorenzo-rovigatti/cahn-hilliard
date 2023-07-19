#include "CahnHilliard.h"
#include "models/Landau.h"

#include <fstream>

int main(int argc, char *argv[]) {
	srand48(51328);

	cxxopts::Options options("cahn-hilliard", "A simple code to simulate spinodal decomposition through the Cahn-Hilliard equation");
	options.add_options()
	("s,steps", "Number of iterations", cxxopts::value<long long int>())
	("l,load-from", "Load the initial conditions from this file", cxxopts::value<std::string>())
	("T,temperature", "Temperature (in Kelvin)", cxxopts::value<double>()->default_value("300"))
	("t,tetramer-density", "Average value of the initial density of tetramers", cxxopts::value<double>()->default_value("5e-5"))
	("noise", "Random noise for the initial psi", cxxopts::value<double>()->default_value("0.01"))
	("R", "Fraction of inter-species linkers", cxxopts::value<double>()->default_value("0.1"))
	("p,print-every", "Number of iterations every which the state of the system will be appended to the trajectory.dat file (0 means never)", cxxopts::value<long long int>()->default_value("0"))
	("h,help", "Print usage");

	auto result = options.parse(argc, argv);

	ch::FreeEnergyModel *model = new ch::Landau(options);
	ch::CahnHilliard<DIM> system(model, options);

	if(argc == 1 || result.count("help")) {
		fprintf(stderr, "%s", options.help().c_str());
		exit(0);
	}

	if(result["steps"].count() == 0) {
		fprintf(stderr, "ERROR: The -s/--steps argument in mandatory\n");
		exit(1);
	}
	long long int steps = result["steps"].as<long long int>();
	long long int print_every = result["print-every"].as<long long int>();

	std::ofstream trajectory[model->N_species()];
	if(print_every > 0) {
		for(int i = 0; i < model->N_species(); i++) {
			char filename[256];
			sprintf(filename, "trajectory_%d.dat", i);
			trajectory[i].open(filename);
		}
	}

	for(int i = 0; i < model->N_species(); i++) {
		char filename[256];
		sprintf(filename, "init_%d.dat", i);
		std::ofstream output(filename);
		system.print_state(i, output);
		output.close();
	}
	system.print_density("initial_density.dat");

	for(long long int t = 0; t < steps; t++) {
		if(print_every > 0 && t % print_every == 0) {
			for(int i = 0; i < model->N_species(); i++) {
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
		for(int i = 0; i < model->N_species(); i++) {
			trajectory[i].close();
		}
	}

	for(int i = 0; i < model->N_species(); i++) {
		char filename[256];
		sprintf(filename, "last_%d.dat", i);
		std::ofstream output(filename);
		system.print_state(i, output);
		output.close();
	}
	system.print_density("final_density.dat");

	return 0;
}
