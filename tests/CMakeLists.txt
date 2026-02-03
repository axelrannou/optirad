# Tests configuration
# find_package(GTest REQUIRED)

# Integration test for L-BFGS optimizer (can run standalone without GTest)
add_executable(test_lbfgs_optimizer
    optimization/test_lbfgs_optimizer.cpp
)

target_include_directories(test_lbfgs_optimizer PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_lbfgs_optimizer PRIVATE optirad_optimization)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_lbfgs_optimizer PRIVATE OpenMP::OpenMP_CXX)
endif()

# Register as CTest
add_test(NAME LBFGSOptimizerTest COMMAND test_lbfgs_optimizer)

# add_subdirectory(core)
# add_subdirectory(geometry)
# add_subdirectory(dose)
# add_subdirectory(optimization)
