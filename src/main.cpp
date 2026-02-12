#include "CahnHilliard.h"
#include "SimulationState.h"
#include "integrators/Integrator.h"
#include "models/Landau.h"
#include "models/GenericWertheim.h"
#include "models/RicciWertheim.h"
#include "models/SalehWertheim.h"
#include "models/SimpleWertheim.h"
#include "utils/Printer.h"
#include "utils/strings.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <memory>
#include <filesystem>

namespace ch {

class Manager final: public Object {
public:
	Manager(toml::table &config) {
		_steps = _config_value<long long int>(config, "steps");

		if(_steps < 0) {
			critical("steps should be a number larger than 0");
		}

		_sim_state.internal_to_user = _config_optional_value<double>(config, "distance_scaling_factor", 1.0);
    	_sim_state.user_to_internal = 1.0 / _sim_state.internal_to_user;

		_print_mass_every = _config_optional_value<long long int>(config, "print_every", 0);

		_print_average_pressure = _config_optional_value<bool>(config, "print_average_pressure", false);
		_print_pressure_every = _config_optional_value<long long int>(config, "print_pressure_every", 0);
		if(_print_pressure_every > 0) {
			_print_average_pressure = true;
		}

		srand48(_config_optional_value<long long int>(config, "seed", std::time(NULL)));

		std::string outp = _config_optional_value<std::string>(config, "output_path", ".");
		_output_path = std::filesystem::path(outp);
		if(!std::filesystem::exists(_output_path) || !std::filesystem::is_directory(_output_path)) {
			critical("Output path '{}' does not exist or is not a directory", outp);
		}

		std::string model_name = _config_value<std::string>(config, "free_energy");
		if(model_name == "landau") {
			_sim_state.model = std::make_unique<ch::Landau>(config);
		}
		else if(model_name == "simple_wertheim") {
			_sim_state.model = std::make_unique<ch::SimpleWertheim>(config);
		}
		else if(model_name == "saleh") {
			_sim_state.model = std::make_unique<ch::SalehWertheim>(config);
		}
		else if(model_name == "generic_wertheim") {
			_sim_state.model = std::make_unique<ch::GenericWertheim>(config);
		}
		else if(model_name == "ricci") {
			_sim_state.model = std::make_unique<ch::RicciWertheim>(config);
		}
		else {
			critical("Unsupported free energy model '{}'", model_name);
		}

		_system = std::make_unique<ch::CahnHilliard<DIM>>(_sim_state, _sim_state.model.get(), config);
		_printer = std::make_unique<ch::Printer<DIM>>(_sim_state, config);

		if(config["load_from"]) {
			_openmode = std::ios_base::app;
			info("Restarting the simulation from a previous configuration: the output information and configurations will be appended to the respective files");
			// we open the file and try to read out the time step 
			std::string filename = _config_value<std::string>(config, "load_from");
			std::ifstream load_from(filename.c_str());

			// the CahnHilliard object has already parsed this file, so we can take for granted
			// that it exists and it is well-formed
			std::string line;
			std::getline(load_from, line);
			utils::trim(line);
			auto spl = utils::split(line);
			if(spl.size() > 3 && spl[0] == "#" && spl[1] == "step" && spl[2] == "=") {
				_initial_t = std::stoll(spl[3]);
				info("Initial time step (as parsed from the 'load_from' file): {}", _initial_t);
			}
			else {
				warning("The header line of the 'load_from' file does not contain information about the time step", _initial_t);
			}

			load_from.close();
		}
	}

	~Manager() {

	}

	void run() {
		_print_current_state("init_", 0);

		std::ofstream mass_output((_output_path / "energy.dat").string(), _openmode);
		std::ofstream pressure_output;
		if(_print_pressure_every > 0) {
			pressure_output.open((_output_path / "pressure.dat").string(), _openmode);
		}

		for(_t = _initial_t; _t < _initial_t + _steps; _t++) {
			if(_printer->should_print_last(_t)) {
				_print_current_state("last_", _t);
			}
			if(_printer->should_print_traj(_t)) {
				for(int i = 0; i < _sim_state.model->N_species(); i++) {
					_sim_state.integrator->sync();
					_printer->add_to_trajectory(i, _t);
				}
				_traj_printed++;
			}
			if(_should_print_pressure(_t)) {
				_system->print_pressure(pressure_output, _t);
			}
			if(_print_mass_every > 0 && _t % _print_mass_every == 0) {
				std::string output_line = fmt::format("{:.5} {:.8} {:.5} {:L}", _t * _system->dt, _system->average_free_energy(), _system->average_mass(), _t);
				if(_print_average_pressure) {
					output_line += fmt::format(" {:.5}", _system->average_pressure());
				}
				mass_output << output_line << std::endl;
				std::cout << output_line << std::endl;
			}
			_system->evolve();
		}

		mass_output.close();
		if(pressure_output.is_open()) {
			pressure_output.close();
		}

		_print_current_state("last_", _steps);
	}

	GET_NAME(Manager)

private:
	void _print_current_state(std::string_view prefix, long long int t) {
		_printer->print_current_state(prefix, t);
		if(_print_average_pressure) {
			_system->print_pressure((_output_path / fmt::format("{}pressure.dat", prefix)).string(), t);
		}
	}

	bool _should_print_pressure(long long int t) {
		return (_print_pressure_every > 0 && t % _print_pressure_every == 0);
	}

	SimulationState<DIM> _sim_state;

	bool _print_average_pressure;
	long long int _initial_t = 0;
	long long int _t;
	long long int _steps, _print_mass_every, _print_pressure_every;
	int _traj_printed = 0;
	std::ios_base::openmode _openmode = std::ios_base::out;

	std::unique_ptr<ch::CahnHilliard<DIM>> _system;
	std::unique_ptr<ch::Printer<DIM>> _printer;

	// directory where output files will be placed (defaults to current folder)
	std::filesystem::path _output_path = std::filesystem::path(".");
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
