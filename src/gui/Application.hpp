#pragma once

#include "Window.hpp"
#include "Renderer.hpp"
#include "panels/IPanel.hpp"
#include <vector>
#include <memory>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

namespace optirad {

class PatientPanel;
class SliceView;
class View3D;

class Application {
public:
    Application();
    ~Application();

    bool init();
    void run();
    void shutdown();

private:
    void renderMenuBar();
    
    GLFWwindow* m_window = nullptr;
    
    // Views
    std::unique_ptr<View3D> m_view3D;
    std::unique_ptr<SliceView> m_axialView;
    std::unique_ptr<SliceView> m_sagittalView;
    std::unique_ptr<SliceView> m_coronalView;
    
    // Panels
    std::unique_ptr<PatientPanel> m_patientPanel;
};

} // namespace optirad
