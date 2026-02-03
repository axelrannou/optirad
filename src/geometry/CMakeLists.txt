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
