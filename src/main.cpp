#include "CahnHilliard.h"
#include "models/Landau.h"

#include <iostream>
#include <fstream>
#include <ctime>

int main(int argc, char *argv[]) {
	if(argc < 2) {
		fprintf(stderr, "Usage is %s configuration_file\n", argv[0]);
		exit(0);
	}

	toml::parse_result res = toml::parse_file(argv[1]);
	if (!res) {
		std::cerr << "Parsing failed with error '" << res.error().description() << "'" << std::endl;
		return 1;
	}

	toml::table config = std::move(res).table();

	if(!config["steps"] || config["steps"].value<long long int>().value() < 0) {
		fprintf(stderr, "steps should be a number larger than 0\n");
		exit(1);
	}
	long long int steps = config["steps"].value<long long int>().value();
	long long int print_every = config["print-every"].value_or(0);
	long long int seed = config["seed"].value_or(0);

	if(config["seed"]) {
		srand48(config["seed"].value<long long int>().value());
	}
	else {
		srand48(std::time(NULL));
	}

	ch::FreeEnergyModel *model = new ch::Landau(config);
	ch::CahnHilliard<DIM> system(model, config);

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
