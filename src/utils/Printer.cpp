/*
 * Printer.cpp
 *
 * Created on: 2/11/2026
 *     Author: Lorenzo
*/

#include "Printer.h"

namespace ch {

template <int dims>
Printer<dims>::Printer(SimulationState<dims> &sim_state, toml::table &config) : _sim_state(sim_state) {
    N = _config_value<int>(config, "N");
	dt = _config_value<double>(config, "dt");
	dx = _config_optional_value<double>(config, "dx", 1.0);

    _print_traj_strategy = _config_optional_value<std::string>(config, "print_trajectory_strategy", "linear");
		
    if(_print_traj_strategy == "linear") {
        _print_trajectory_every = _config_optional_value<long long int>(config, "print_trajectory_every", 0);
        _print_last_every = _config_optional_value<long long int>(config, "print_last_every", _print_trajectory_every);
    }
    else if(_print_traj_strategy == "log") {
        _log_n0 = _config_value<int>(config, "log_n0");
        _log_fact = _config_value<double>(config, "log_fact");
        _print_last_every = _config_value<long long int>(config, "print_last_every");
    }
    else {
        critical("Unsupported printing strategy '{}'", _print_traj_strategy);
    }

    std::string outp = _config_optional_value<std::string>(config, "output.path", ".");
    _output_path = std::filesystem::path(outp);
    if(!std::filesystem::exists(_output_path)) {
        critical("Output path '{}' does not exist", outp);
    }
    if(!std::filesystem::is_directory(_output_path)) {
        critical("Output path '{}' is not a directory", outp);
    }

    std::string trajp;
    _print_vtk = _config_optional_value<bool>(config, "output.print_vtk", false);
    // trajectory_path is mandatory if we print vtk trajectories, otherwise it defaults to output.path
    if(_print_vtk) {
        trajp = _config_value<std::string>(config, "output.trajectory_path");
        info("VTK trajectory printing is enabled, trajectory files will be written to '{}'", trajp);
    }
    else {
        trajp = _config_optional_value<std::string>(config, "output.trajectory_path", outp);
        info("Trajectory files will be written in native format to '{}'", trajp);
    }

    _traj_path = std::filesystem::path(trajp);
    if(!std::filesystem::exists(_traj_path)) {
        critical("Trajectory output path '{}' does not exist", trajp);
    }
    if(!std::filesystem::is_directory(_traj_path)) {
        critical("Trajectory output path '{}' is not a directory", trajp);
    }

    // if we print native trajectories and we are not loading from a previous state, then we
    // open the trajectory files in output mode (overwriting any existing file with the same name)
    if(!_print_vtk && !config["load_from"]) {
        for(int i = 0; i < _sim_state.model->N_species(); i++) {
            std::string filename = (_traj_path / fmt::format("traj_{}.dat", i)).string();
            std::ofstream output(filename, std::ios_base::out);
            output.close();
        }
    }
}

template <int dims>
Printer<dims>::~Printer() {

}

template <int dims>
void Printer<dims>::print_current_state(std::string_view prefix, long long int time_step) {
    for(int i = 0; i < _sim_state.model->N_species(); i++) {
        auto filename = (_output_path / fmt::format("{}{}.dat", prefix, i)).string();
        _write_native(filename, i, time_step);

        if(_print_vtk) {
            auto vtk_filename = (_output_path / fmt::format("{}{}.vtk", prefix, i)).string();
            _write_vtk(vtk_filename, i, time_step);
        }
    }

    auto filename = (_output_path / fmt::format("{}density.dat", prefix)).string();
    _write_native(filename, -1, time_step); // -1 indicates that we want to print the total density

    if(_print_vtk) {
        auto vtk_filename = (_output_path / fmt::format("{}density.vtk", prefix)).string();
        _write_vtk(vtk_filename, -1, time_step);
    }
}

template <int dims>
void Printer<dims>::add_to_trajectory(int species, long long int time_step) {
    if(_print_vtk) {
        auto vtk_filename = (_traj_path / fmt::format("traj_{}_{}.vtk", species, time_step)).string();
        _write_vtk(vtk_filename, species, time_step);
    }
    else {
        std::string filename = (_traj_path / fmt::format("traj_{}.dat", species)).string();
        std::ofstream output(filename, std::ios_base::app);
        _write_native(output, species, time_step);
        output.close();
    }
}

template <int dims>
bool Printer<dims>::should_print_last(long long int time_step) {
    return (_print_last_every > 0 && time_step % _print_last_every == 0);
}

template <int dims>
bool Printer<dims>::should_print_traj(long long int time_step) {
    if(_print_traj_strategy == "linear") {
        return (_print_trajectory_every > 0 && time_step % _print_trajectory_every == 0);
    }
    else if(_print_traj_strategy == "log") {
        long long int next_t = (long long int) round((_log_n0 * std::pow(_log_fact, _traj_printed)));
        return (next_t == time_step);
    }
    return false;
}

template <int dims>
void Printer<dims>::_write_native(std::ofstream &output, int species, long long int time_step) {
    const int nx = N;
    const int ny = (dims >= 2) ? N : 1;
    const int nz = (dims >= 3) ? N : 1;
    const int grid_size = nx * ny * nz;

    output << fmt::format("# step = {}, t = {:.5}, size = {}x{}x{}", time_step, time_step * dt, nx, ny, nz) << std::endl;
	int newline_every = ny;
	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0 && idx % newline_every == 0) {
			output << std::endl;
		}
        double rho_value = (species == -1) ? _sim_state.rho.field_sum(idx) : _sim_state.rho(idx, species);
		output << _sim_state.density_to_user(rho_value) << " ";
	}
	output << std::endl;
}

template <int dims>
void Printer<dims>::_write_native(const std::string &filename, int species, long long int time_step) {
    std::ofstream output(filename);
    _write_native(output, species, time_step);
    output.close();
}

template <int dims>
void Printer<dims>::_write_vtk(const std::string &filename, int species, long long int time_step) {
    std::ofstream output(filename);

    const int nx = N;
    const int ny = (dims >= 2) ? N : 1;
    const int nz = (dims >= 3) ? N : 1;
    const int grid_size = nx * ny * nz;
    const double sx = _sim_state.length_to_user(dx);
    const double sy = _sim_state.length_to_user(dx);
    const double sz = _sim_state.length_to_user(dx);

    output << "# vtk DataFile Version 3.0\n";
    output << "CIRCA scalar output\n";
    output << "ASCII\n";
    output << "DATASET STRUCTURED_POINTS\n";
    output << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    output << "ORIGIN 0 0 0\n";
    output << "SPACING " << std::setprecision(16) << sx << " " << sy << " " << sz << "\n";
    output << "POINT_DATA " << (static_cast<size_t>(nx) * ny * nz) << "\n";
    output << "SCALARS species_" << species << " double 1\n";
    output << "LOOKUP_TABLE default\n";

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                int idx = i;
                if constexpr (dims >= 2) {
                    idx += j * nx;
                }
                if constexpr (dims >= 3) {
                    idx += k * nx * ny;
                }
                double rho_value = (species == -1) ? _sim_state.rho.field_sum(idx) : _sim_state.rho(idx, species);
                output << std::setprecision(16) << _sim_state.density_to_user(rho_value) << "\n";
            }
        }
    }

    output.close();
}

template class Printer<1>;
template class Printer<2>;
template class Printer<3>;

} /* namespace ch */
