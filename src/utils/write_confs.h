#pragma once
#include <fstream>
#include <iomanip>
#include <string>

namespace ch::utils {

namespace detail {
template <int D>
int dim_or_one(int grid_size, int idx) {
    return (idx < D) ? grid_size : 1;
}
template <int D>
double dx_or_one(double dx, int idx) {
    return (idx < D) ? dx : 1.0;
}
}  // namespace detail

template <int dims>
void write_ch(SimulationState &sim_state, std::ofstream &output, int species, int grid_size, double dx, long long int time_step, double dt) {
    const int nx = detail::dim_or_one<dims>(grid_size, 0);
    const int ny = detail::dim_or_one<dims>(grid_size, 1);
    const int nz = detail::dim_or_one<dims>(grid_size, 2);
    const double sx = sim_state.length_to_user(detail::dx_or_one<dims>(dx, 0));
    const double sy = sim_state.length_to_user(detail::dx_or_one<dims>(dx, 1));
    const double sz = sim_state.length_to_user(detail::dx_or_one<dims>(dx, 2));

    output << fmt::format("# step = {}, t = {:.5}, size = {}x{}x{}", time_step, time_step * dt, nx, ny, nz) << std::endl;
	int newline_every = (dims == 1) ? 1 : grid_size;
	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0 && idx % newline_every == 0) {
			output << std::endl;
		}
		output << sim_state.density_to_user(sim_state.rho(idx, species)) << " ";
	}
	output << std::endl;
}

template <int D>
inline int flat(const std::array<int, D>& I, const int N) {
    int s = 0, m = 1;
    for (int d = 0; d < D; ++d) {
        s += I[d] * m;
        m *= N;
    }
    return s;
}

// Write a single scalar field to VTK (STRUCTURED_POINTS, ASCII).
// Works for D=1/2/3; for D<3 we set missing dims to 1 and nz=1 (2D) or ny=nz=1 (1D).
template <int D>
void write_vtk(SimulationState &sim_state, std::ofstream &output, int species, int N, double dx, long long int time_step, double dt) {
    const int nx = detail::dim_or_one<D>(N, 0);
    const int ny = detail::dim_or_one<D>(N, 1);
    const int nz = detail::dim_or_one<D>(N, 2);
    const double sx = sim_state.length_to_user(detail::dx_or_one<D>(dx, 0));
    const double sy = sim_state.length_to_user(detail::dx_or_one<D>(dx, 1));
    const double sz = sim_state.length_to_user(detail::dx_or_one<D>(dx, 2));

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

    // VTK expects x fastest, then y, then z. Our flat() does x-fastest too,
    // so we can index with flat({i,j,k}).
    for(int k = 0; k < nz; ++k) {
        for(int j = 0; j < ny; ++j) {
            for(int i = 0; i < nx; ++i) {
                std::array<int, 3> I3{i, j, k};
                // Build an index array of size D: missing dims are 0.
                std::array<int, D> ID{};
                if constexpr (D >= 1) {
                    ID[0] = I3[0];
                }
                if constexpr (D >= 2) {
                    ID[1] = I3[1];
                }
                if constexpr (D >= 3) {
                    ID[2] = I3[2];
                }
                const int idx = flat<D>(ID, N);
                output << std::setprecision(16) << sim_state.density_to_user(sim_state.rho(idx, species)) << "\n";
            }
        }
    }
}

}  // namespace ch::utils
