/*
 * defs_CUDA.h
 *
 *  Created on: Nov 19, 2023
 *      Author: lorenzo
 */

#ifndef SRC_DEFS_CUDA_H_
#define SRC_DEFS_CUDA_H_

#ifndef NOCUDA

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <stdio.h>

/// CUDA_SAFE_CALL replacement for backwards compatibility (CUDA < 5.0)
#define CUDA_SAFE_CALL(call)                                  \
  do {                                                        \
    cudaError_t err = call;                                   \
    if (err != cudaSuccess) {                                 \
      printf("CUDA error at %s %d: %s\n", __FILE__, __LINE__, \
             cudaGetErrorString(err));                        \
      exit(EXIT_FAILURE);                                     \
    }                                                         \
  } while (0)
/// CUT_CHECK_ERROR replacement for backwards compatibility (CUDA < 5.0)
#define CUT_CHECK_ERROR(x) getLastCudaError(x);

#ifndef CUFFT_CALL
#define CUFFT_CALL( call )                                                                                             \
    {                                                                                                                  \
        auto status = static_cast<cufftResult>( call );                                                                \
        if ( status != CUFFT_SUCCESS )                                                                                 \
            fprintf( stderr,                                                                                           \
                     "ERROR: CUFFT call \"%s\" in line %d of file %s failed "                                          \
                     "with "                                                                                           \
                     "code (%d).\n",                                                                                   \
                     #call,                                                                                            \
                     __LINE__,                                                                                         \
                     __FILE__,                                                                                         \
                     status );                                                                                         \
    }
#endif

#define BLOCK_SIZE 64

/// threads per block
#define TINBLOCK (blockDim.x*blockDim.y)
/// c_number of blocks
#define NBLOCKS (gridDim.x*gridDim.y)
/// c_number of threads
#define NTHREADS (NBLOCKS * TINBLOCK)

/// thread id relative to its block
#define TID (blockDim.x*threadIdx.y + threadIdx.x)
/// block id
#define BID (gridDim.x*blockIdx.y + blockIdx.x)
/// thread id
#define IND (TINBLOCK * BID + TID)

#define COPY_ARRAY_TO_CONSTANT(dest, src, size) {\
		float *val = new float[(size)];\
		for(int i = 0; i < (size); i++) val[i] = (float) ((src)[i]);\
		CUDA_SAFE_CALL(cudaMemcpyToSymbol((dest), val, (size)*sizeof(float)));\
		delete[] val; }

#define COPY_NUMBER_TO_FLOAT(dest, src) {\
		float tmp = src;\
		CUDA_SAFE_CALL(cudaMemcpyToSymbol((dest), &tmp, sizeof(float)));\
		}

#ifdef CUDA_FIELD_FLOAT
using cufftFieldComplex = cufftComplex;
#else
using cufftFieldComplex = cufftDoubleComplex;
#endif /* CUDA_FIELD_FLOAT */

#endif /* NOCUDA */

#ifdef CUDA_FIELD_FLOAT
using field_type = float;
#else
using field_type = double;
#endif /* CUDA_FIELD_FLOAT */

#endif /* SRC_DEFS_CUDA_H_ */