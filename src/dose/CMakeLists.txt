add_library(optirad_dose
    DoseMatrix.cpp
    DoseInfluenceMatrix.cpp
    DoseEngineFactory.cpp
    engines/PencilBeamEngine.cpp
)

target_include_directories(optirad_dose PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_dose PUBLIC optirad_geometry)
