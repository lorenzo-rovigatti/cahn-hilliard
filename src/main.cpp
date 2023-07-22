#include "CahnHilliard.h"
#include "models/Landau.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <memory>

namespace ch {

class Manager final: public Object {
public:
	Manager(toml::table &config) {
		if(!config["steps"] || config["steps"].value<long long int>().value() < 0) {
			critical("steps should be a number larger than 0");
		}

		_steps = config["steps"].value<long long int>().value();
		_print_mass_every = config["print-every"].value<long long int>().value_or(0);
		_print_trajectory_every = config["print-trajectory-every"].value<long long int>().value_or(0);

		if(config["seed"]) {
			srand48(config["seed"].value<long long int>().value());
		}
		else {
			srand48(std::time(NULL));
		}

		_model = std::make_unique<ch::Landau>(config);
		_system = std::make_unique<ch::CahnHilliard<DIM>>(_model.get(), config);

		_trajectories.resize(_model->N_species());
		if(_print_trajectory_every > 0) {
			for(int i = 0; i < _model->N_species(); i++) {
				_trajectories[i].open(fmt::format("trajectory_{}.dat", i));
			}
		}
	}

	~Manager() {
		for(auto &traj : _trajectories) {
			if(traj.good()) {
				traj.close();
			}
		}
	}

	void run() {
		_print_current_state("init_");

		std::ofstream mass_output("mass.dat");

		for(long long int t = 0; t < _steps; t++) {
			if(_print_trajectory_every > 0 && t % _print_trajectory_every == 0) {
				for(int i = 0; i < _model->N_species(); i++) {
					_system->print_species_density(i, _trajectories[i]);
					_system->print_species_density(i, fmt::format("last_{}.dat", i));
				}
			}
			if(_print_mass_every > 0 && t % _print_mass_every == 0) {
				std::string mass_line = fmt::format("{} {} {}", t * _system->dt, _system->total_mass(), t);
				mass_output << mass_line << std::endl;
				std::cout << mass_line << std::endl;
			}
			_system->evolve();
		}

		mass_output.close();

		_print_current_state("last_");
	}

	GET_NAME(Manager)

private:
	void _print_current_state(std::string_view prefix) {
		for(int i = 0; i < _model->N_species(); i++) {
			_system->print_species_density(i, fmt::format("{}{}.dat", prefix, i));
		}
		_system->print_total_density(fmt::format("{}density.dat", prefix));
	}

	long long int _steps, _print_mass_every, _print_trajectory_every;

	std::vector<std::ofstream> _trajectories;
	std::unique_ptr<ch::FreeEnergyModel> _model;
	std::unique_ptr<ch::CahnHilliard<DIM>> _system;
};

}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		std::cerr << fmt::format("Usage is {} configuration_file", argv[0]) << std::endl;
		return 0;
	}

	toml::parse_result res = toml::parse_file(argv[1]);
	if(!res) {
		std::cerr << fmt::format("Parsing failed with error '{}'", res.error().description()) << std::endl;
		return 1;
	}

	toml::table config = std::move(res).table();

	ch::Manager manager(config);
	manager.run();

	return 0;
}
