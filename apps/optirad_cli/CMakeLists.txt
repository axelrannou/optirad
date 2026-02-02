add_executable(optirad_cli main.cpp)

target_link_libraries(optirad_cli PRIVATE
    optirad_optimization
    optirad_io
)
