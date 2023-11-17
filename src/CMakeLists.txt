set(common_files
	CahnHilliard.cpp
	Object.cpp
	models/Landau.cpp
	models/SalehWertheim.cpp
	models/SimpleWertheim.cpp
	utils/Delta.cpp
	)

IF(CUDA)
	CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
	FIND_PACKAGE("CUDA" REQUIRED)

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

	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -use_fast_math --compiler-options "-ffast-math -O3" -Xcudafe "--diag_suppress=unrecognized_gcc_pragma")

	LIST(APPEND common_SOURCES
		CUDA/CUDA_sort.cu
		CUDA/CUDAUtils.cu
		CUDA/cuda_utils/cuda_device_utils.cu
		CUDA/Backends/CUDABaseBackend.cu
		CUDA/Backends/MD_CUDABackend.cu
		CUDA/Interactions/CUDABaseInteraction.cu
		CUDA/Interactions/CUDADNAInteraction.cu
		CUDA/Interactions/CUDALJInteraction.cu
		CUDA/Interactions/CUDARNAInteraction.cu
		CUDA/Interactions/CUDAPatchyInteraction.cu
		CUDA/Interactions/CUDATEPInteraction.cu
		CUDA/Interactions/CUDAInteractionFactory.cu
		CUDA/Lists/CUDAListFactory.cu
		CUDA/Lists/CUDANoList.cu
		CUDA/Lists/CUDASimpleVerletList.cu
		CUDA/Lists/CUDABinVerletList.cu
		CUDA/Thermostats/CUDABaseThermostat.cu
		CUDA/Thermostats/CUDANoThermostat.cu
		CUDA/Thermostats/CUDAThermostatFactory.cu
		CUDA/Thermostats/CUDABrownianThermostat.cu
		CUDA/Thermostats/CUDASRDThermostat.cu
		CUDA/Thermostats/CUDALangevinThermostat.cu
		CUDA/Thermostats/CUDABussiThermostat.cu
		Utilities/Timings.cpp
	)

	IF(CUDA_DOUBLE)
		ADD_DEFINITIONS(-DCUDA_DOUBLE_PRECISION)
	ELSE()
		LIST(APPEND common_SOURCES
			CUDA/Backends/MD_CUDAMixedBackend.cu 
			CUDA/Backends/FFS_MD_CUDAMixedBackend.cu
		)
	ENDIF()

	# we need the library to be shared or the oxpy compilation fails
	CUDA_ADD_LIBRARY(${lib_name} SHARED ${common_SOURCES})
	CUDA_ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES})
	CUDA_ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES} Utilities/Timings.cpp)

	SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
endif(CUDA)
	
add_library(ch_common ${common_files})
	
add_executable(ch_1D main.cpp ${source_files})
target_link_libraries(ch_1D ch_common spdlog::spdlog)
target_compile_definitions(ch_1D PRIVATE DIM=1)

add_executable(ch_2D main.cpp ${source_files})
target_link_libraries(ch_2D ch_common spdlog::spdlog)
target_compile_definitions(ch_2D PRIVATE DIM=2)

add_compile_options(-Wall -Wextra)
