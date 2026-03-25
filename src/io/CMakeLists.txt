add_library(optirad_io
    DicomImporter.cpp
    RTStructParser.cpp
    MachineLoader.cpp
)

target_include_directories(optirad_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_io PUBLIC optirad_core optirad_geometry optirad_segmentation optirad_dose nlohmann_json::nlohmann_json)
target_compile_definitions(optirad_io PRIVATE OPTIRAD_DATA_DIR="${OPTIRAD_DATA_DIR}")

# DCMTK for DICOM handling (optional - stub implementation if not found)
find_package(DCMTK QUIET)
if(DCMTK_FOUND)
    target_link_libraries(optirad_io PRIVATE ${DCMTK_LIBRARIES})
    target_include_directories(optirad_io PRIVATE ${DCMTK_INCLUDE_DIRS})
    target_compile_definitions(optirad_io PRIVATE OPTIRAD_HAS_DCMTK)
endif()

# TBB for parallel structure parsing (optional - falls back to sequential)
find_package(TBB QUIET)
if(TBB_FOUND)
    target_link_libraries(optirad_io PRIVATE TBB::tbb)
    target_compile_definitions(optirad_io PRIVATE OPTIRAD_HAS_TBB)
    message(STATUS "TBB found - parallel RT-STRUCT parsing enabled")
else()
    message(STATUS "TBB not found - RT-STRUCT parsing will use sequential mode")
endif()
