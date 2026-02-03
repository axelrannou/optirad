add_library(optirad_optimization
    ObjectiveFunction.cpp
    Constraint.cpp
    OptimizerFactory.cpp
    objectives/SquaredDeviation.cpp
    objectives/SquaredOverdose.cpp
    objectives/SquaredUnderdose.cpp
    objectives/DVHObjective.cpp
    optimizers/LBFGSOptimizer.cpp
)

target_include_directories(optirad_optimization PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(optirad_optimization PUBLIC optirad_dose optirad_geometry optirad_utils)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_optimization PUBLIC OpenMP::OpenMP_CXX)
endif()
