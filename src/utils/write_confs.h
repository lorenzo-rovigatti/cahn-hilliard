#pragma once
#include <fstream>
#include <iomanip>
#include <string>

namespace ch::utils {

namespace detail {
template <int D>
int dim_or_one(const std::array<int, D>& n, int idx) {
    return (idx < D) ? n[idx] : 1;
}
template <int D>
double dx_or_one(const std::array<double, D>& dx, int idx) {
    return (idx < D) ? dx[idx] : 1.0;
}
}  // namespace detail

template <int dims>
void write_ch(SimulationState &sim_state, std::ofstream &output, int species, int grid_size, double dx) {
    const int nx = detail::dim_or_one<dims>(grid_size, 0);
    const int ny = detail::dim_or_one<dims>(grid_size, 1);
    const int nz = detail::dim_or_one<dims>(grid_size, 2);
    const double sx = detail::dx_or_one<dims>(dx, 0);
    const double sy = detail::dx_or_one<dims>(dx, 1);
    const double sz = detail::dx_or_one<dims>(dx, 2);

    double dt = 1; // CHANGE THIS: we should pass the actual dt here, but for now we just want to print the time in physical units, so we can use dt=1 and multiply by it when printing the time

    output << fmt::format("# step = {}, t = {:.5}, size = {}x{}x{}x", sim_state.time_step, sim_state.time_step * dt, nx, ny, nz) << std::endl;
	int newline_every = (dims == 1) ? 1 : N;
	for(int idx = 0; idx < grid_size; idx++) {
		if(idx > 0 && idx % newline_every == 0) {
			output << std::endl;
		}
        double density = _sim_state.rho(idx, species) / (sim_state.internal_to_user * sim_state.internal_to_user * sim_state.internal_to_user);
		output << _density_to_user() << " ";
	}
	output << std::endl;
}

template <int D>
inline int flat(const std::array<int, D>& I, const int grid_size) {
    int s = 0, m = 1;
    for (int d = 0; d < D; ++d) {
        s += I[d] * m;
        m *= grid_size;
    }
    return s;
}

// Write a single scalar field to VTK (STRUCTURED_POINTS, ASCII).
// Works for D=1/2/3; for D<3 we set missing dims to 1 and nz=1 (2D) or ny=nz=1 (1D).
template <int D>
void write_vtk(SimulationState &state, std::ofstream &output, int species, int grid_size, double dx) {
    const int nx = detail::dim_or_one<D>(grid_size, 0);
    const int ny = detail::dim_or_one<D>(grid_size, 1);
    const int nz = detail::dim_or_one<D>(grid_size, 2);
    const double sx = detail::dx_or_one<D>(dx, 0);
    const double sy = detail::dx_or_one<D>(dx, 1);
    const double sz = detail::dx_or_one<D>(dx, 2);

    output << "# vtk DataFile Version 3.0\n";
    output << "CIRCA scalar output\n";
    output << "ASCII\n";
    output << "DATASET STRUCTURED_POINTS\n";
    output << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    output << "ORIGIN 0 0 0\n";
    output << "SPACING " << std::setprecision(16) << sx << " " << sy << " " << sz << "\n";
    output << "POINT_DATA " << (static_cast<size_t>(nx) * ny * nz) << "\n";
    output << "SCALARS " << scalar_name << " double 1\n";
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
                const int lin = flat<D>(ID, grid_size);
                output << std::setprecision(16) << f.a[lin] << "\n";
            }
        }
    }
}

}  // namespace ch::utils
