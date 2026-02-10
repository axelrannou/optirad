add_library(optirad_io
    DicomImporter.cpp
    RTStructParser.cpp
    DicomExporter.cpp
)

target_include_directories(optirad_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_io PUBLIC optirad_core optirad_geometry)

# DCMTK for DICOM handling (optional - stub implementation if not found)
find_package(DCMTK QUIET)
if(DCMTK_FOUND)
    target_link_libraries(optirad_io PRIVATE ${DCMTK_LIBRARIES})
    target_include_directories(optirad_io PRIVATE ${DCMTK_INCLUDE_DIRS})
    target_compile_definitions(optirad_io PRIVATE OPTIRAD_HAS_DCMTK)
endif()
