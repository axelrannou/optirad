add_library(optirad_geometry
    CoordinateSystem.cpp
    Grid.cpp
    Volume.cpp
    Structure.cpp
    StructureSet.cpp
)

target_include_directories(optirad_geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_geometry PUBLIC optirad_core)
