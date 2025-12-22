set(common_files
	CahnHilliard.cpp
	Object.cpp
	integrators/BailoFiniteVolume.cpp
	integrators/EulerCPU.cpp
	integrators/EulerMobilityCPU.cpp
	integrators/GelMobilityCPU.cpp
	integrators/Integrator.cpp
	integrators/PseudospectralCPU.cpp
	mobility/IMobility.cpp
	models/Landau.cpp
	models/GenericWertheim.cpp
	models/RicciWertheim.cpp
	models/SalehWertheim.cpp
	models/SimpleWertheim.cpp
	utils/Delta.cpp
	utils/strings.cpp
	../extern/fsolve/fsolve.cpp
	)

IF(CUDA)
	CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
	FIND_PACKAGE("CUDA" REQUIRED)

	SET(CMAKE_CUDA_STANDARD 17)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

	INCLUDE_DIRECTORIES(${CUDA_INSTALL_PREFIX}/include/)
	LINK_DIRECTORIES(${CUDA_INSTALL_PREFIX}/lib)

	INCLUDE(../select_compute_arch.cmake)

	IF(CUDA_COMMON_ARCH)
		CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Common)
	ELSE(CUDA_COMMON_ARCH)
		CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
	ENDIF(CUDA_COMMON_ARCH)

	MESSAGE(STATUS "Generating code for the following CUDA architectures: ${ARCH_FLAGS_readable}")
	LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
	MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

	IF(G)
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lineinfo)
	ENDIF(G)

	IF(DEBUG)
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lineinfo -G)
	ENDIF(DEBUG)

	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; --expt-relaxed-constexpr -use_fast_math --compiler-options "-ffast-math -O3 --std=c++17" -Xcudafe "--diag_suppress=unrecognized_gcc_pragma")

	LIST(APPEND common_files
		CUDA/integrators/CUDAIntegrator.cu
		CUDA/integrators/EulerCUDA.cu
		CUDA/integrators/EulerMobilityCUDA.cu
		CUDA/integrators/PseudospectralCUDA.cu
		CUDA/models/GenericWertheim.cu
		CUDA/models/Landau.cu
		CUDA/models/RicciWertheim.cu
		CUDA/models/SalehWertheim.cu
		CUDA/models/SimpleWertheim.cu
	)

	cuda_add_library(ch_common ${common_files})
	set_property(TARGET ch_common PROPERTY CUDA_STANDARD 17)

	find_package(CUDAToolkit REQUIRED)
	target_link_libraries(ch_common CUDA::cufft)

	cuda_add_executable(ch_1D main.cpp)
	cuda_add_executable(ch_2D main.cpp)

	set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
else()
	add_definitions(-DNOCUDA)

	add_library(ch_common ${common_files})
	add_executable(ch_1D main.cpp)
	add_executable(ch_2D main.cpp)
endif(CUDA)
	
target_link_libraries(ch_1D ch_common spdlog::spdlog fftw3)
target_compile_definitions(ch_1D PRIVATE DIM=1)

target_link_libraries(ch_2D ch_common spdlog::spdlog fftw3)
target_compile_definitions(ch_2D PRIVATE DIM=2)

add_compile_options(-Wall -Wextra)
