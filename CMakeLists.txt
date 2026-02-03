cmake_minimum_required(VERSION 3.20)
project(optirad VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(OPTIRAD_BUILD_TESTS "Build unit tests" ON)
option(OPTIRAD_BUILD_GUI "Build GUI application" ON)

# External dependencies
add_subdirectory(external)

# Source libraries
add_subdirectory(src)

# Applications
add_subdirectory(apps)

# Tests
if(OPTIRAD_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
