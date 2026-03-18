# Check if OpenGL and GLFW are available
find_package(OpenGL)
find_package(glfw3)
find_package(GLEW REQUIRED)

if(NOT OpenGL_FOUND OR NOT glfw3_FOUND)
    message(WARNING "OpenGL or GLFW not found. GUI will not be built.")
    message(STATUS "To install on Ubuntu/Debian: sudo apt-get install libgl1-mesa-dev libglfw3-dev")
    message(STATUS "To install on Fedora/RHEL: sudo dnf install mesa-libGL-devel glfw-devel")
    return()
endif()

# Check if Dear ImGui is available
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
if(NOT EXISTS ${IMGUI_DIR}/imgui.h)
    message(WARNING "Dear ImGui not found at ${IMGUI_DIR}")
    message(STATUS "To download: cd ${CMAKE_SOURCE_DIR}/external && git clone --branch docking --single-branch https://github.com/ocornut/imgui.git")
    return()
endif()

# Dear ImGui library
add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)
target_link_libraries(imgui PUBLIC OpenGL::GL glfw)

# GUI library
add_library(optirad_gui
    Application.cpp
    Theme.cpp
    stb_image_impl.cpp
    Window.cpp
    Renderer.cpp
    panels/PatientPanel.cpp
    panels/BeamPanel.cpp
    panels/OptimizationPanel.cpp
    panels/DoseStatsPanel.cpp
    panels/LogPanel.cpp
    panels/PlanningPanel.cpp
    panels/StfPanel.cpp
    panels/PhaseSpacePanel.cpp
    views/SliceView.cpp
    views/DVHView.cpp
    views/View3D.cpp
    views/renderers/CubeRenderer.cpp
    views/renderers/AxisLabels.cpp
    views/renderers/VolumeRenderer.cpp
    views/renderers/StructureRenderer.cpp
    views/renderers/BeamRenderer.cpp
    views/renderers/PhaseSpaceRenderer.cpp
)

target_include_directories(optirad_gui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(optirad_gui PRIVATE ${CMAKE_SOURCE_DIR}/external)

target_link_libraries(optirad_gui PUBLIC
    optirad_core
    optirad_io
    optirad_geometry
    optirad_optimization
    optirad_phsp
    imgui
    OpenGL::GL
    glfw
    GLEW::GLEW
)

# Make OPTIRAD_DATA_DIR available for file browser default path
target_compile_definitions(optirad_gui PRIVATE OPTIRAD_DATA_DIR="${OPTIRAD_DATA_DIR}")
# Path to gui image assets (theme icon, logos)
target_compile_definitions(optirad_gui PRIVATE OPTIRAD_GUI_IMG_DIR="${CMAKE_CURRENT_SOURCE_DIR}/img")

# OpenMP (needed for header-only steering module used by StfPanel)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(optirad_gui PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "GUI: OpenMP found, STF generation will be parallelised")
endif()

message(STATUS "GUI library will be built with Dear ImGui")
