/*
 * EulerMobilityCUDA.h
 *
 * Created on: 3/29/2024
 *      Author: Lorenzo
*/

#ifndef EULERMOBILITYCUDA_H
#define EULERMOBILITYCUDA_H

#include "../../defs_CUDA.h"

#include <curand_kernel.h>

#include "CUDAIntegrator.h"

namespace ch {

template <int dims>
struct CUDAVector {
    std::array<field_type, dims> values;

    __host__ __device__ field_type& operator[](int i) { return values[i]; }
    __host__ __device__ const field_type& operator[](int i) const { return values[i]; }

    __host__ __device__ CUDAVector<dims> operator+(const CUDAVector<dims>& other) const {
        CUDAVector<dims> result;
        for (int i = 0; i < dims; i++) {
            result[i] = values[i] + other[i];
        }
        return result;
    }

    __host__ __device__ CUDAVector<dims> operator-(const CUDAVector<dims>& other) const {
        CUDAVector<dims> result;
        for (int i = 0; i < dims; i++) {
            result[i] = values[i] - other[i];
        }
        return result;
    }

    __host__ __device__ CUDAVector<dims>& operator+=(const CUDAVector<dims>& other) {
        for (int i = 0; i < dims; i++) {
            values[i] += other[i];
        }
        return *this;
    }

    __host__ __device__ CUDAVector<dims>& operator-=(const CUDAVector<dims>& other) {
        for (int i = 0; i < dims; i++) {
            values[i] -= other[i];
        }
        return *this;
    }
};

template <int dims, typename DataType>
class CUDAGrid {
public:
    std::array<int, dims> sizes;
    int bins;
    int species;
    int species_size;
    int total_size;
    DataType *d_data;

    CUDAGrid() {

    }

    CUDAGrid(int mbins, int mspecies) : bins(mbins), species(mspecies) {
        species_size = 1;
        for(int i = 0; i < dims; i++) {
            sizes[i] = bins;
            species_size *= sizes[i];
        }
        total_size = species_size * species;
        CUDA_SAFE_CALL(cudaMalloc(&d_data, total_size * sizeof(DataType)));
    }

    // Destructor
    ~CUDAGrid() {
        CUDA_SAFE_CALL(cudaFree(d_data));
    }

    // Convert N-dimensional indices to 1D index (with periodic boundary conditions)
    __host__ __device__ int index(const std::array<int, dims> &indices, int species) const {
        int idx = 0, stride = 1;
        for (int i = 0; i < dims; i++) {
            idx += ((indices[i] + sizes[i]) % sizes[i]) * stride;
            stride *= sizes[i];
        }
        idx += species_size * species;
        return idx;
    }

    // Get device pointer
    __host__ __device__ DataType* data() { return d_data; }
};

template<int dims>
class EulerMobilityCUDA : public CUDAIntegrator<dims> {
public:
    EulerMobilityCUDA(SimulationState<dims> &sim_state,FreeEnergyModel *model, toml::table &config);

    ~EulerMobilityCUDA();

    void evolve() override;

    GET_NAME(EulerMobilityCUDA)

protected:
    bool _with_noise = false;
    CUDAGrid<dims, CUDAVector<dims>> *_h_flux = nullptr;
    CUDAGrid<dims, CUDAVector<dims>> *_d_flux = nullptr;
    curandState *_d_rand_states = nullptr;

    bool _supports_nonconstant_mobility() const override {
        return true;
    }
};

} /* namespace ch */

#endif /* EULERMOBILITYCUDA_H */
