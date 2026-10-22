add_library(optirad_dose
    DoseMatrix.cpp
    DoseInfluenceMatrix.cpp
    DoseEngineFactory.cpp
    DijSerializer.cpp
    SiddonRayTracer.cpp
    SSDCalculator.cpp
    RadDepthCalculator.cpp
    FluenceDoseModel.cpp
    FluenceDoseOperator.cpp
    engines/PencilBeamEngine.cpp
)

target_include_directories(optirad_dose PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(optirad_dose PUBLIC optirad_geometry optirad_core)
target_compile_definitions(optirad_dose PRIVATE OPTIRAD_DATA_DIR="${OPTIRAD_DATA_DIR}")

option(OPTIRAD_ENABLE_CUDA "Build CUDA matrix-free forward and adjoint dose backend" OFF)
if(OPTIRAD_ENABLE_CUDA)
    if(CMAKE_VERSION VERSION_LESS 3.18)
        message(FATAL_ERROR "OPTIRAD_ENABLE_CUDA requires CMake 3.18 or newer")
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_sources(optirad_dose PRIVATE CudaFluenceBackend.cu CudaRayTracer.cu CudaDirectDose.cu)
    target_link_libraries(optirad_dose PRIVATE CUDA::cudart CUDA::cufft)
    target_compile_definitions(optirad_dose PUBLIC OPTIRAD_HAS_CUDA)
    set_target_properties(optirad_dose PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set_property(TARGET optirad_dose PROPERTY CUDA_ARCHITECTURES 75)
    endif()
    set_source_files_properties(CudaFluenceBackend.cu PROPERTIES COMPILE_OPTIONS "--fmad=false")
    set_source_files_properties(CudaRayTracer.cu PROPERTIES COMPILE_OPTIONS "--fmad=false")
    set_source_files_properties(CudaDirectDose.cu PROPERTIES COMPILE_OPTIONS "--fmad=false")
endif()

# OpenMP for parallel dose computation (optional)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_dose PRIVATE OpenMP::OpenMP_CXX)
endif()
