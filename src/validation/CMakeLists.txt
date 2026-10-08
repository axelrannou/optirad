add_library(optirad_validation
    GammaAnalysis.cpp
    DoseRawIO.cpp
    WaterPhantom.cpp
)

target_include_directories(optirad_validation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(optirad_validation PUBLIC optirad_workflow optirad_io optirad_dose optirad_core optirad_geometry)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_validation PRIVATE OpenMP::OpenMP_CXX)
endif()
