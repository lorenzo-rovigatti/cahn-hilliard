/*
 * Printer.h
 *
 * Created on: 2/11/2026
 *      Author: Lorenzo
*/

#ifndef PRINTER_H
#define PRINTER_H

#include "../Object.h"
#include "../SimulationState.h"

#include <filesystem>
#include <string>

namespace ch {

class Printer : public Object {
public:
    Printer(SimulationState &sim_state, toml::table &config);
    ~Printer();

    template <int dims>
    void print_current_state(std::string_view prefix, long long int t);

    GET_NAME(Printer)

private:
    bool _should_print_last(long long int t);
    bool _should_print_traj(long long int t);

    template <int dims>
    void _write_native(std::ofstream &output, int species, long long int time_step);

    template <int dims>
    void _write_native(const std::string &filename, int species, long long int time_step);

    template <int dims>
    void _write_vtk(const std::string &filename, int species, long long int time_step);

    SimulationState &_sim_state;

    int N, grid_size;
    double dt;
    double dx;

    std::string _print_traj_strategy;
	long long int _print_trajectory_every, _print_last_every;
    bool _print_vtk;
	int _traj_printed = 0;
	int _log_n0 = 0;
	double _log_fact = 0;
	std::ios_base::openmode _openmode = std::ios_base::out;
    // directory where output files will be placed (defaults to current folder)
	std::filesystem::path _output_path;
    // directory where trajectory files will be placed (defaults to _output_path)
	std::filesystem::path _traj_path;
};

} /* namespace ch */

#endif /* PRINTER_H */
