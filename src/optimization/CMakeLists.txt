add_library(optirad_optimization
    ObjectiveFunction.cpp
    Constraint.cpp
    OptimizerFactory.cpp
    objectives/SquaredDeviation.cpp
    objectives/SquaredOverdose.cpp
    objectives/SquaredUnderdose.cpp
    optimizers/LBFGSOptimizer.cpp
)

target_include_directories(optirad_optimization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_optimization PUBLIC optirad_dose)
