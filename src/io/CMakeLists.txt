add_library(optirad_io
    DicomImporter.cpp
    DicomExporter.cpp
)

target_include_directories(optirad_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_io PUBLIC optirad_geometry)
# target_link_libraries(optirad_io PRIVATE dcmtk::dcmtk)
