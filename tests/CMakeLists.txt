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

# Machine loader test
add_executable(test_machine_loader
    io/test_machine_loader.cpp
)
target_include_directories(test_machine_loader PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_machine_loader PRIVATE optirad_io GTest::gtest)
target_compile_definitions(test_machine_loader PRIVATE OPTIRAD_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_test(NAME MachineLoaderTest COMMAND test_machine_loader)

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
target_link_libraries(test_photon_imrt_stf_generator PRIVATE optirad_io optirad_core optirad_geometry GTest::gtest GTest::gtest_main)
target_compile_definitions(test_photon_imrt_stf_generator PRIVATE OPTIRAD_DATA_DIR="${OPTIRAD_DATA_DIR}")
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

# ============================================================================
# UNIT TESTS (GoogleTest) - Phase-Space Module
# ============================================================================

add_executable(test_iaea_header_parser
    phsp/test_iaea_header_parser.cpp
)
target_include_directories(test_iaea_header_parser PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_iaea_header_parser PRIVATE optirad_phsp GTest::gtest GTest::gtest_main)
target_compile_definitions(test_iaea_header_parser PRIVATE OPTIRAD_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_test(NAME IAEAHeaderParserTest COMMAND test_iaea_header_parser)

add_executable(test_iaea_phsp_reader
    phsp/test_iaea_phsp_reader.cpp
)
target_include_directories(test_iaea_phsp_reader PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_iaea_phsp_reader PRIVATE optirad_phsp GTest::gtest GTest::gtest_main)
target_compile_definitions(test_iaea_phsp_reader PRIVATE OPTIRAD_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_test(NAME IAEAPhspReaderTest COMMAND test_iaea_phsp_reader)

add_executable(test_phase_space_data
    phsp/test_phase_space_data.cpp
)
target_include_directories(test_phase_space_data PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_phase_space_data PRIVATE optirad_phsp GTest::gtest GTest::gtest_main)
add_test(NAME PhaseSpaceDataTest COMMAND test_phase_space_data)

add_executable(test_phase_space_beam_source
    phsp/test_phase_space_beam_source.cpp
)
target_include_directories(test_phase_space_beam_source PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_phase_space_beam_source PRIVATE optirad_phsp optirad_core GTest::gtest GTest::gtest_main)
target_compile_definitions(test_phase_space_beam_source PRIVATE OPTIRAD_DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_test(NAME PhaseSpaceBeamSourceTest COMMAND test_phase_space_beam_source)

# ============================================================================
# UNIT TESTS (GoogleTest) - Dose Module
# ============================================================================

add_executable(test_dose_influence_matrix
    dose/test_dose_influence_matrix.cpp
)
target_include_directories(test_dose_influence_matrix PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_dose_influence_matrix PRIVATE optirad_dose GTest::gtest GTest::gtest_main)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_dose_influence_matrix PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME DoseInfluenceMatrixTest COMMAND test_dose_influence_matrix)

add_executable(test_dij_serializer
    dose/test_dij_serializer.cpp
)
target_include_directories(test_dij_serializer PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_dij_serializer PRIVATE optirad_dose GTest::gtest GTest::gtest_main)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_dij_serializer PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME DijSerializerTest COMMAND test_dij_serializer)

add_executable(test_fft2d
    dose/test_fft2d.cpp
)
target_include_directories(test_fft2d PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_fft2d PRIVATE optirad_dose GTest::gtest GTest::gtest_main)
add_test(NAME FFT2DTest COMMAND test_fft2d)

add_executable(test_grid_interpolant2d
    dose/test_grid_interpolant2d.cpp
)
target_include_directories(test_grid_interpolant2d PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_grid_interpolant2d PRIVATE optirad_dose GTest::gtest GTest::gtest_main)
add_test(NAME GridInterpolant2DTest COMMAND test_grid_interpolant2d)

add_executable(test_siddon_ray_tracer
    dose/test_siddon_ray_tracer.cpp
)
target_include_directories(test_siddon_ray_tracer PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_siddon_ray_tracer PRIVATE optirad_dose optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME SiddonRayTracerTest COMMAND test_siddon_ray_tracer)

add_executable(test_ssd_calculator
    dose/test_ssd_calculator.cpp
)
target_include_directories(test_ssd_calculator PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_ssd_calculator PRIVATE optirad_dose optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME SSDCalculatorTest COMMAND test_ssd_calculator)

add_executable(test_rad_depth_calculator
    dose/test_rad_depth_calculator.cpp
)
target_include_directories(test_rad_depth_calculator PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_rad_depth_calculator PRIVATE optirad_dose optirad_geometry GTest::gtest GTest::gtest_main)
add_test(NAME RadDepthCalculatorTest COMMAND test_rad_depth_calculator)

add_executable(test_dose_engine
    dose/test_dose_engine.cpp
)
target_include_directories(test_dose_engine PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_dose_engine PRIVATE optirad_dose optirad_geometry optirad_core GTest::gtest GTest::gtest_main)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_dose_engine PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME DoseEngineTest COMMAND test_dose_engine)
