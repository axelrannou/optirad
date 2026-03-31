#pragma once

#include "Window.hpp"
#include "Renderer.hpp"
#include "AppState.hpp"
#include "Theme.hpp"
#include "panels/IPanel.hpp"
#include <vector>
#include <memory>
#include <string>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

namespace optirad {

class PatientPanel;
class PlanningPanel;
class StfPanel;
class PhaseSpacePanel;
class OptimizationPanel;
class DoseStatsPanel;
class DVHPanel;
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
    void render3DViewWindow();
    void setupDockLayout();
    /// Load a PNG file into an OpenGL RGBA texture. Returns false on failure.
    static bool loadTexture(const std::string& path, GLuint* outTexture);
    
    GLFWwindow* m_window = nullptr;
    
    // DPI scale factor
    float m_dpiScale = 1.0f;

    // Theme
    AppTheme m_theme = AppTheme::Dark;
    GLuint m_themeIconTexture = 0;

    // Layout state
    bool m_layoutInitialized = false;
    bool m_view3DVisible = true;
    
    // Shared state
    GuiAppState m_appState;

    // Views
    std::unique_ptr<View3D> m_view3D;
    std::unique_ptr<SliceView> m_axialView;
    std::unique_ptr<SliceView> m_sagittalView;
    std::unique_ptr<SliceView> m_coronalView;
    
    // Panels
    std::unique_ptr<PatientPanel> m_patientPanel;
    std::unique_ptr<PlanningPanel> m_planningPanel;
    std::unique_ptr<StfPanel> m_stfPanel;
    std::unique_ptr<PhaseSpacePanel> m_phaseSpacePanel;
    std::unique_ptr<OptimizationPanel> m_optimizationPanel;
    std::unique_ptr<DoseStatsPanel> m_doseStatsPanel;
    std::unique_ptr<DVHPanel> m_dvhPanel;

    // Track whether imported dose was already added to DoseManager
    bool m_importedDoseAdded = false;

    bool focus3DView = false;
};

} // namespace optirad
