cmake_minimum_required(VERSION 3.15)
project(optirad CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Data directory (machine files, etc.)
set(OPTIRAD_DATA_DIR "${CMAKE_SOURCE_DIR}/data" CACHE PATH "Path to OptiRad data directory")

# Options
option(OPTIRAD_BUILD_TESTS "Build unit tests" ON)
option(OPTIRAD_BUILD_GUI "Build GUI application" ON)

# External dependencies
add_subdirectory(external)

# Core modules (always built)
add_subdirectory(src/utils)
add_subdirectory(src/core)
add_subdirectory(src/geometry)
add_subdirectory(src/io)
add_subdirectory(src/dose)
add_subdirectory(src/optimization)
add_subdirectory(src/phsp)

# GUI (optional - only if dependencies found)
add_subdirectory(src/gui)

# Applications
add_subdirectory(apps)

# Tests
if(OPTIRAD_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
