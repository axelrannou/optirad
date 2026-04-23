add_library(optirad_dose
    DoseMatrix.cpp
    DoseInfluenceMatrix.cpp
    DoseEngineFactory.cpp
    DijSerializer.cpp
    SiddonRayTracer.cpp
    SSDCalculator.cpp
    RadDepthCalculator.cpp
    engines/PencilBeamEngine.cpp
)

target_include_directories(optirad_dose PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(optirad_dose PUBLIC optirad_geometry optirad_core)
target_compile_definitions(optirad_dose PRIVATE OPTIRAD_DATA_DIR="${OPTIRAD_DATA_DIR}")

# OpenMP for parallel dose computation (optional)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_dose PRIVATE OpenMP::OpenMP_CXX)
endif()
