add_library(optirad_geometry
    Grid.cpp
    Volume.cpp
    Structure.cpp
    StructureSet.cpp
    CoordinateSystem.cpp
)

target_include_directories(optirad_geometry PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(optirad_geometry PUBLIC optirad_core)

# OpenMP for parallel voxel dilation (optional - falls back to sequential)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_geometry PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found - parallel voxel dilation enabled")
else()
    message(STATUS "OpenMP not found - voxel dilation will use sequential mode")
endif()

