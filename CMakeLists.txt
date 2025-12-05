cmake_minimum_required(VERSION 3.10)

option(CUDA "Set to ON to compile with CUDA support" OFF)
option(Debug "Set to ON to compile with debug symbols" OFF)
option(G "Set to ON to compile with optimisations and debug symbols" OFF)
option(NATIVE_COMPILATION "Set to OFF to compile without the -march=native flag. This may be required when compiling binaries to be used elsewhere" ON)
option(CUDA_COMMON_ARCH "Set to OFF to autodetect the GPU and compile for its architecture, set to ON (default) to compile for the most common architectures" ON)

project(ch)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

if(Debug)
	set(CMAKE_BUILD_TYPE Debug)
	message(STATUS "Compiling in Debug mode")
endif(Debug)

if(G)
	set(CMAKE_BUILD_TYPE RelWithDebInfo)
	message(STATUS "Compiling in Release+Debug symbols mode")
endif(G)

if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
	add_definitions(-ffast-math)
	
	if(NATIVE_COMPILATION)
		include(CheckCXXCompilerFlag)
        unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
        CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
                ADD_DEFINITIONS(-march=native)
		else()
			message(STATUS "The '-march=native' flag is not supported, compiling without it")
        endif()
	else()
		message(STATUS "Compiling without the '-march=native' flag")
	endif()
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# build spdlog as a shared lib
add_subdirectory(extern/spdlog)
set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

include_directories(extern extern/spdlog/include)
add_subdirectory(src)

# add tests
include(CTest)
enable_testing()

add_test(
    NAME python_tests
    COMMAND pytest -q
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
