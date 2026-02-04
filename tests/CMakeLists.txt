# Tests configuration

# Integration test for L-BFGS optimizer
add_executable(test_lbfgs_optimizer
    optimization/test_lbfgs_optimizer.cpp
)
target_include_directories(test_lbfgs_optimizer PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_lbfgs_optimizer PRIVATE optirad_optimization)

# DICOM importer test
add_executable(test_dicom_importer
    io/test_dicom_importer.cpp
)
target_include_directories(test_dicom_importer PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_dicom_importer PRIVATE optirad_io)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_lbfgs_optimizer PRIVATE OpenMP::OpenMP_CXX)
endif()

# Register tests
add_test(NAME LBFGSOptimizerTest COMMAND test_lbfgs_optimizer)
add_test(NAME DicomImporterTest COMMAND test_dicom_importer ${CMAKE_SOURCE_DIR}/data/dicom)
