#include "../defs_CUDA.h"

template<int dims> 
__device__ void _fill_coords(int N, int bits, int coords[dims], int idx) {
    for(int d = 0; d < dims; d++) {
		coords[d] = idx & (N - 1);
		idx >>= bits; // divide by N
	}
}

template<int dims>
__device__ int _cell_idx(int bits, int coords[dims]) {
	int idx = 0;
	int multiply_by = 1;
	for(int d = 0; d < dims; d++) {
		idx += coords[d] * multiply_by;
		multiply_by <<= bits; // multiply by N
	}
	return idx;
}

template<int dims, typename number>
__device__ float _cell_laplacian_old(int size, int N, int bits, number *field, int idx, float dx) {
    int base_idx = (idx / size) * size;
    int rel_idx = idx % size;
    int N_minus_one = N - 1;

    float laplacian = 0.f;

    if constexpr (dims == 1) {
        int rel_idx_m = (rel_idx - 1 + N) & N_minus_one;
        int rel_idx_p = (rel_idx + 1) & N_minus_one;

        laplacian = ((float) field[base_idx + rel_idx_m] + (float) field[base_idx + rel_idx_p] - 2.f * (float) field[idx]) / (dx * dx);
    }
    else if constexpr (dims == 2) {
        int coords_xy[2];
        _fill_coords<2>(N, bits, coords_xy, rel_idx);

        int coords_xmy[2] = {
                (coords_xy[0] - 1 + N) & N_minus_one,
                coords_xy[1]
        };

        int coords_xym[2] = {
                coords_xy[0],
                (coords_xy[1] - 1 + N) & N_minus_one
        };

        int coords_xpy[2] = {
                (coords_xy[0] + 1) & N_minus_one,
                coords_xy[1]
        };

        int coords_xyp[2] = {
                coords_xy[0],
                (coords_xy[1] + 1) & N_minus_one
        };

        laplacian = (
                (float) field[base_idx + _cell_idx<2>(bits, coords_xmy)] +
                (float) field[base_idx + _cell_idx<2>(bits, coords_xpy)] +
                (float) field[base_idx + _cell_idx<2>(bits, coords_xym)] +
                (float) field[base_idx + _cell_idx<2>(bits, coords_xyp)] -
                4.f * (float) field[idx])
                / (dx * dx);
    }
    else if constexpr (dims == 3) {
        int coords_xyz[3];
        _fill_coords<3>(N, bits, coords_xyz, rel_idx);

        int coords_xmyz[3] = {
                (coords_xyz[0] - 1 + N) & N_minus_one,
                coords_xyz[1],
                coords_xyz[2]
        };

        int coords_xyzm[3] = {
                coords_xyz[0],
                (coords_xyz[1] - 1 + N) & N_minus_one,
                coords_xyz[2]
        };

        int coords_xyzm2[3] = {
                coords_xyz[0],
                coords_xyz[1],
                (coords_xyz[2] - 1 + N) & N_minus_one
        };

        int coords_xpyz[3] = {
                (coords_xyz[0] + 1) & N_minus_one,
                coords_xyz[1],
                coords_xyz[2]
        };

        int coords_xyzp[3] = {
                coords_xyz[0],
                (coords_xyz[1] + 1) & N_minus_one,
                coords_xyz[2]
        };

        int coords_xyzp2[3] = {
                coords_xyz[0],
                coords_xyz[1],
                (coords_xyz[2] + 1) & N_minus_one
        };

        laplacian = (
                (float) field[base_idx + _cell_idx<3>(bits, coords_xmyz)] +
                (float) field[base_idx + _cell_idx<3>(bits, coords_xpyz)] +
                (float) field[base_idx + _cell_idx<3>(bits, coords_xyzm)] +
                (float) field[base_idx + _cell_idx<3>(bits, coords_xyzp)] +
                (float) field[base_idx + _cell_idx<3>(bits, coords_xyzm2)] +
                (float) field[base_idx + _cell_idx<3>(bits, coords_xyzp2)] -
                6.f * (float) field[idx])
                / (dx * dx);
    }
    else {
        static_assert(dims == 1 || dims == 2 || dims == 3, "Unsupported dims for the laplacian");
    }

    return laplacian;
}

template<int dims, typename number>
__device__ float _cell_laplacian(int size, int N, int bits, number *field, int idx, float dx) {
    int base_idx = (idx / size) * size;
    int rel_idx = idx % size;
    int N_minus_one = N - 1;

    int coords[dims];
    _fill_coords<dims>(N, bits, coords, rel_idx);

    float laplacian = 0.f;

    // Iterate through each dimension to compute neighbor contributions
    #pragma unroll
    for(int d = 0; d < dims; d++) {
        int coords_p[dims], coords_m[dims];

        #pragma unroll
        for(int i = 0; i < dims; i++) {
            coords_p[i] = coords[i];
            coords_m[i] = coords[i];
        }

        // Positive and negative neighbors in dimension d
        coords_p[d] = (coords[d] + 1) & N_minus_one;
        coords_m[d] = (coords[d] - 1 + N) & N_minus_one;

        laplacian += (float) field[base_idx + _cell_idx<dims>(bits, coords_p)];
        laplacian += (float) field[base_idx + _cell_idx<dims>(bits, coords_m)];
    }

    // Subtract center value: (2 * dims) times
    laplacian -= 2.f * dims * (float) field[idx];
    laplacian /= (dx * dx);

    return laplacian;
}