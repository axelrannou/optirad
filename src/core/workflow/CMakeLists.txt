add_library(optirad_workflow
    PlanBuilder.cpp
    DoseCalculationPipeline.cpp
    PhaseSpaceBuilder.cpp
    OptimizationPipeline.cpp
)

target_include_directories(optirad_workflow PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(optirad_workflow PUBLIC
    optirad_core
    optirad_io
    optirad_dose
    optirad_optimization
    optirad_phsp
    optirad_geometry
    optirad_utils
)

# OpenMP for PhaseSpaceBuilder parallel beam loading + steering (header-only)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_workflow PRIVATE OpenMP::OpenMP_CXX)
endif()
