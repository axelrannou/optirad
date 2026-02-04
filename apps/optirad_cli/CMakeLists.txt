add_executable(optirad_cli main.cpp)

target_include_directories(optirad_cli PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(optirad_cli PRIVATE 
    optirad_io 
    optirad_core 
    optirad_geometry
    optirad_utils
)
