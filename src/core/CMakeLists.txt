add_library(optirad_core
    Patient.cpp
    Plan.cpp
    Beam.cpp
    Machine.cpp
)

target_include_directories(optirad_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_core PUBLIC optirad_utils)
