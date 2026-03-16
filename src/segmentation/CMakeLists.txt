add_library(optirad_segmentation
    BodyContourGenerator.cpp
)

target_include_directories(optirad_segmentation PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(optirad_segmentation PUBLIC optirad_geometry)

# TBB for parallel per-slice processing (optional - falls back to sequential)
find_package(TBB QUIET)
if(TBB_FOUND)
    target_link_libraries(optirad_segmentation PRIVATE TBB::tbb)
    target_compile_definitions(optirad_segmentation PRIVATE OPTIRAD_HAS_TBB)
    message(STATUS "TBB found - parallel body contour generation enabled")
else()
    message(STATUS "TBB not found - body contour generation will use sequential mode")
endif()
