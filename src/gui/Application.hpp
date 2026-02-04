#pragma once

#include "Window.hpp"
#include "Renderer.hpp"
#include "panels/IPanel.hpp"
#include <vector>
#include <memory>
#include <GLFW/glfw3.h>

namespace optirad {

class PatientPanel;
class SliceView;

class Application {
public:
    Application();
    ~Application();

    bool init();
    void run();
    void shutdown();

private:
    void update();
    void render();
    void renderMenuBar();

    GLFWwindow* m_window = nullptr;  // Raw pointer for GLFW

    // Panels
    std::unique_ptr<PatientPanel> m_patientPanel;

    // Views
    std::unique_ptr<SliceView> m_axialView;
    std::unique_ptr<SliceView> m_sagittalView;
    std::unique_ptr<SliceView> m_coronalView;
};

} // namespace optirad
