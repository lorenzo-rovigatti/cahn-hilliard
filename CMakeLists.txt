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

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# build spdlog as a shared lib
add_subdirectory(extern/spdlog)
set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

include_directories(extern extern/spdlog/include)
add_subdirectory(src)
