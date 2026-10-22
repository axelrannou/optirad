add_library(optirad_sequencing
    LeafSequencer.cpp
    DeliverableDoseCalculator.cpp
    ApertureFluence.cpp
    ArcSequencer.cpp
    SlidingWindow.cpp
)

target_include_directories(optirad_sequencing PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(optirad_sequencing PUBLIC optirad_core optirad_dose optirad_geometry)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_sequencing PRIVATE OpenMP::OpenMP_CXX)
endif()
