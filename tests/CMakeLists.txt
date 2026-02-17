# Tests configuration
enable_testing()

# Find GoogleTest
find_package(GTest REQUIRED)

# ============================================================================
# LEGACY TESTS (Manual main() functions) - For reference/compatibility
# ============================================================================

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

# ============================================================================
# UNIT TESTS (GoogleTest) - Core Module
# ============================================================================

add_executable(test_patient
    core/test_patient.cpp
)
target_include_directories(test_patient PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_patient PRIVATE optirad_core GTest::gtest GTest::gtest_main)
add_test(NAME PatientTest COMMAND test_patient)

add_executable(test_plan
    core/test_plan.cpp
)
target_include_directories(test_plan PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_plan PRIVATE optirad_core GTest::gtest GTest::gtest_main)
add_test(NAME PlanTest COMMAND test_plan)

add_executable(test_beam
    core/test_beam.cpp
)
target_include_directories(test_beam PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_beam PRIVATE optirad_core GTest::gtest GTest::gtest_main)
add_test(NAME BeamTest COMMAND test_beam)

# Ray test
add_executable(test_ray
    core/test_ray.cpp
)
target_include_directories(test_ray PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_ray PRIVATE optirad_core GTest::gtest GTest::gtest_main)
add_test(NAME RayTest COMMAND test_ray)

# ============================================================================
# UNIT TESTS (GoogleTest) - Steering Module
# ============================================================================

add_executable(test_photon_imrt_stf_generator
    steering/test_photon_imrt_stf_generator.cpp
)
target_include_directories(test_photon_imrt_stf_generator PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_photon_imrt_stf_generator PRIVATE optirad_core optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME PhotonIMRTStfGeneratorTest COMMAND test_photon_imrt_stf_generator)

# ============================================================================
# UNIT TESTS (GoogleTest) - Geometry Module
# ============================================================================

add_executable(test_grid
    geometry/test_grid.cpp
)
target_include_directories(test_grid PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_grid PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME GridTest COMMAND test_grid)

add_executable(test_volume
    geometry/test_volume.cpp
)
target_include_directories(test_volume PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_volume PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME VolumeTest COMMAND test_volume)

add_executable(test_coordinate_system
    geometry/test_coordinate_system.cpp
)
target_include_directories(test_coordinate_system PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_coordinate_system PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME CoordinateSystemTest COMMAND test_coordinate_system)

add_executable(test_indexing
    geometry/test_indexing.cpp
)
target_include_directories(test_indexing PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_indexing PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME IndexingTest COMMAND test_indexing)

add_executable(test_voxel_dilation
    geometry/test_voxel_dilation.cpp
)
target_include_directories(test_voxel_dilation PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_voxel_dilation PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME VoxelDilationTest COMMAND test_voxel_dilation)

# ============================================================================
# UNIT TESTS (GoogleTest) - Utils Module
# ============================================================================

add_executable(test_logger
    utils/test_logger.cpp
)
target_include_directories(test_logger PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_logger PRIVATE optirad_core GTest::gtest GTest::gtest_main)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_logger PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME LoggerTest COMMAND test_logger)

add_executable(test_mathutils
    utils/test_mathutils.cpp
)
target_include_directories(test_mathutils PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_mathutils PRIVATE optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME MathUtilsTest COMMAND test_mathutils)

# ============================================================================
# LEGACY TEST REGISTRATION (Manual tests)
# ============================================================================

add_test(NAME LBFGSOptimizerTest COMMAND test_lbfgs_optimizer)
add_test(NAME DicomImporterTest COMMAND test_dicom_importer ${CMAKE_SOURCE_DIR}/data/dicom)
