#include "CahnHilliard.h"
#include "models/Landau.h"

#include <iostream>
#include <fstream>
#include <ctime>

namespace ch {

class Manager final: public Object {
public:
	Manager(toml::table &config) {
		if(!config["steps"] || config["steps"].value<long long int>().value() < 0) {
			critical("steps should be a number larger than 0");
		}

		_steps = config["steps"].value<long long int>().value();
		_print_every = config["print-every"].value<long long int>().value_or(0);

		if(config["seed"]) {
			srand48(config["seed"].value<long long int>().value());
		}
		else {
			srand48(std::time(NULL));
		}

		_model = new ch::Landau(config);
		_system = new ch::CahnHilliard<DIM>(_model, config);

		_trajectories.resize(_model->N_species());
		if(_print_every > 0) {
			for(int i = 0; i < _model->N_species(); i++) {
				char filename[256];
				sprintf(filename, "trajectory_%d.dat", i);
				_trajectories[i].open(filename);
			}
		}
	}

	~Manager() {
		for(auto &traj : _trajectories) {
			if(traj.good()) {
				traj.close();
			}
		}

		delete _system;
		delete _model;
	}

	void run() {
		_print_current_state("init_");

		for(long long int t = 0; t < _steps; t++) {
			if(_print_every > 0 && t % _print_every == 0) {
				for(int i = 0; i < _model->N_species(); i++) {
					_system->print_state(i, _trajectories[i]);

					char filename[256];
					sprintf(filename, "last_%d.dat", i);
					std::ofstream output(filename);
					_system->print_state(i, output);
					output.close();
				}

				fprintf(stdout, "%lld %lf %lf\n", t, t * _system->dt, _system->total_mass());
			}
			_system->evolve();
		}

		_print_current_state("last_");
	}

	GET_NAME(Manager)

private:
	void _print_current_state(std::string_view prefix) {
		for(int i = 0; i < _model->N_species(); i++) {
			std::ofstream output(fmt::format("{}{}.dat", prefix, i));
			_system->print_state(i, output);
			output.close();
		}
		_system->print_density(fmt::format("{}density.dat", prefix));
	}

	long long int _steps, _print_every;

	std::vector<std::ofstream> _trajectories;
	ch::FreeEnergyModel *_model;
	ch::CahnHilliard<DIM> *_system;
};

}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		fprintf(stderr, "Usage is %s configuration_file\n", argv[0]);
		exit(0);
	}

	toml::parse_result res = toml::parse_file(argv[1]);
	if(!res) {
		std::cerr << "Parsing failed with error '" << res.error().description() << "'" << std::endl;
		return 1;
	}

	toml::table config = std::move(res).table();

	ch::Manager manager(config);
	manager.run();

	return 0;
}
