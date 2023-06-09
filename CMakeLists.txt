cmake_minimum_required(VERSION 3.10)

option(Debug "Set to ON to compile with debug symbols" OFF)
option(G "Set to ON to compile with optimisations and debug symbols" OFF)

project(ch)

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
endif()

include_directories(extern)
add_executable(ch ch.cpp)
add_executable(ch_ldm ch_ldm.cpp)
add_executable(ch_binary_ldm ch_ldm_binary.cpp)

add_compile_options(-Wall -Wextra)

add_executable(ch_fft ch_fft.cpp)
add_executable(ch_ldm_binary_fft ch_ldm_binary_fft.cpp)

target_link_libraries(ch_fft fftw3)
target_link_libraries(ch_ldm_binary_fft fftw3)
