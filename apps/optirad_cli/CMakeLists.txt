add_executable(optirad_cli main.cpp)

target_include_directories(optirad_cli PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(optirad_cli PRIVATE 
    optirad_workflow
    optirad_io 
    optirad_core 
    optirad_geometry
    optirad_dose
    optirad_optimization
    optirad_phsp
    optirad_utils
)

# OpenMP for parallel STF generation in header-only steering module
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_cli PRIVATE OpenMP::OpenMP_CXX)
endif()

