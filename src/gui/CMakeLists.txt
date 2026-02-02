add_library(optirad_gui
    Application.cpp
    Window.cpp
    Renderer.cpp
    panels/PatientPanel.cpp
    panels/BeamPanel.cpp
    panels/OptimizationPanel.cpp
    panels/DoseStatsPanel.cpp
    panels/LogPanel.cpp
    views/SliceView.cpp
    views/DVHView.cpp
)

target_include_directories(optirad_gui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(optirad_gui PUBLIC optirad_optimization optirad_io)
# target_link_libraries(optirad_gui PRIVATE imgui glfw glad)
