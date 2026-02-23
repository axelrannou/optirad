# Phase-space module: IAEA file parsing, particle data model, beam source
add_library(optirad_phsp
    IAEAHeaderParser.cpp
    IAEAPhspReader.cpp
    PhaseSpaceData.cpp
    PhaseSpaceBeamSource.cpp
)

target_include_directories(optirad_phsp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

# OpenMP for parallel multi-file reading in PhaseSpaceBeamSource::build()
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    message(STATUS "PHSP: OpenMP found — parallel PSF file reading enabled")
endif()

target_link_libraries(optirad_phsp PUBLIC
    optirad_core
    optirad_geometry
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)
