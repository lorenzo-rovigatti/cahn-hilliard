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
#include <set>

namespace ch {

template<int dims>
class Printer : public Object {
public:
    Printer(SimulationState<dims> &sim_state, toml::table &config);
    ~Printer();

    bool should_print_last(long long int t);
    bool should_print_traj(long long int t);
    bool should_print_pressure(long long int t);

    void write_native(std::ofstream &output, 
        const std::string &obs_name, 
        const MultiField<double> &field, 
        double converting_factor, 
        int species, 
        long long int time_step);

    void write_native(const std::string &filename, 
        const std::string &obs_name, 
        const MultiField<double> &field, 
        double converting_factor, 
        int species, 
        long long int time_step);

    void write_vtk(const std::string &filename, 
        const std::string &obs_name, 
        const MultiField<double> &field, 
        double converting_factor, 
        int species, 
        long long int time_step);

    void print_current_state(std::string_view prefix, long long int time_step);

    void add_to_trajectory(const std::string &prefix,
        const std::string &obs_name,
        const MultiField<double> &field,
        double converting_factor, 
        int species,
        long long int time_step);

    GET_NAME(Printer)

private:
    SimulationState<dims> &_sim_state;

    int N, grid_size;
    double dt;
    double dx;

    std::set<std::string> _valid_trajectory_prefixes;

    std::string _print_traj_strategy;
    int _traj_printed = 0;
	long long int _print_trajectory_every;

    bool _print_pressure;
    std::string _print_pressure_strategy;
	long long int _print_pressure_every;

    long long int _print_last_every;

    bool _print_vtk;
	int _log_n0 = 0;
	double _log_fact = 0;
	std::ios_base::openmode _openmode = std::ios_base::out;
    // directory where output files will be placed (defaults to current folder)
	std::filesystem::path _output_path;
    // directory where trajectory files will be placed (defaults to _output_path)
	std::filesystem::path _traj_path;

    void _check_prefix(const std::string &prefix);
};

} /* namespace ch */

#endif /* PRINTER_H */
