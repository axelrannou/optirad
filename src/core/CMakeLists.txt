add_library(optirad_core
    Patient.cpp
    Plan.cpp
    Beam.cpp
    Machine.cpp
    PatientData.cpp
    Stf.cpp
)

target_include_directories(optirad_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_core PUBLIC optirad_utils optirad_geometry)

# OpenMP for parallel ray generation (optional - falls back to sequential)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_core PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found - parallel STF generation enabled")
else()
    message(STATUS "OpenMP not found - STF generation will use sequential mode")
endif()

add_subdirectory(workflow)
