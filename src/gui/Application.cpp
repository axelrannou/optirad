#include "Application.hpp"
#include "panels/PatientPanel.hpp"
#include "panels/PlanningPanel.hpp"
#include "panels/StfPanel.hpp"
#include "views/SliceView.hpp"
#include "views/View3D.hpp"
#include "utils/Logger.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace optirad {

Application::Application() = default;
Application::~Application() = default;

bool Application::init() {
    // Initialize GLFW
    if (!glfwInit()) {
        Logger::error("Failed to initialize GLFW");
        return false;
    }
    
    // Create window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA
    m_window = glfwCreateWindow(1920, 1080, "OptiRad TPS", nullptr, nullptr);
    if (!m_window) {
        Logger::error("Failed to create GLFW window");
        glfwTerminate();
        return false;
    }
    
    glfwMakeContextCurrent(m_window);
    glfwSwapInterval(1); // VSync

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        Logger::error("Failed to initialize GLEW");
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }

    glEnable(GL_MULTISAMPLE); // Enable MSAA
    
    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    
    if (!ImGui_ImplGlfw_InitForOpenGL(m_window, true)) {
        Logger::error("Failed to initialize ImGui GLFW backend");
        ImGui::DestroyContext();
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }
    
    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        Logger::error("Failed to initialize ImGui OpenGL3 backend");
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }
    
    // Create 3D view
    m_view3D = std::make_unique<View3D>();
    m_view3D->init();
    
    // Set up scroll callback for 3D view zoom
    glfwSetWindowUserPointer(m_window, m_view3D.get());
    glfwSetScrollCallback(m_window, [](GLFWwindow* win, double xOffset, double yOffset) {
        auto* view = static_cast<View3D*>(glfwGetWindowUserPointer(win));
        if (view && !ImGui::GetIO().WantCaptureMouse) {
            view->handleScroll(yOffset);
        }
    });
    
    // Create panels
    m_patientPanel = std::make_unique<PatientPanel>();
    m_planningPanel = std::make_unique<PlanningPanel>(m_appState);
    m_stfPanel = std::make_unique<StfPanel>(m_appState);
    
    // Create slice views
    m_axialView = std::make_unique<SliceView>(SliceOrientation::Axial);
    m_sagittalView = std::make_unique<SliceView>(SliceOrientation::Sagittal);
    m_coronalView = std::make_unique<SliceView>(SliceOrientation::Coronal);
    
    return true;
}

void Application::run() {
    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();
        
        m_view3D->handleMouseInput(m_window);
        
        // Clear buffers
        int display_w, display_h;
        glfwGetFramebufferSize(m_window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Update 3D view with patient data
        if (m_patientPanel->hasData()) {
            auto* data = m_patientPanel->getPatientData();

            // Bridge PatientPanel data → GuiAppState (shared_ptr from raw ptr)
            if (!m_appState.dicomLoaded()) {
                // Wrap PatientPanel's data in a non-owning shared_ptr (PatientPanel still owns it)
                m_appState.patientData = std::shared_ptr<PatientData>(
                    data, [](PatientData*) {}); // no-op deleter
            }

            m_view3D->setPatientData(data);
            m_axialView->setPatientData(data);
            m_sagittalView->setPatientData(data);
            m_coronalView->setPatientData(data);
        }

        // Pass STF to 3D view when available
        if (m_appState.stfGenerated()) {
            m_view3D->setStf(m_appState.stf.get());
            m_stfPanel->setBeamRenderer(m_view3D->getBeamRenderer());
        }
        
        // Render 3D view
        m_view3D->render();
        
        // Start ImGui frame (renders on top of 3D scene)
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // Render main menu
        renderMenuBar();
        
        // Render panels
        m_patientPanel->render();
        m_planningPanel->render();
        m_stfPanel->render();
        
        // Render views
        m_axialView->render();
        m_sagittalView->render();
        m_coronalView->render();
        
        // Render ImGui
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        glfwSwapBuffers(m_window);
    }
}

void Application::shutdown() {
    // Cleanup 3D view
    if (m_view3D) {
        m_view3D->cleanup();
    }
    
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    
    if (m_window) {
        glfwDestroyWindow(m_window);
    }
    glfwTerminate();
}

void Application::renderMenuBar() {
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Import DICOM...")) {
                // Trigger import dialog in patient panel
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) {
                glfwSetWindowShouldClose(m_window, true);
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("View")) {
            ImGui::MenuItem("Patient Panel", nullptr, true);
            ImGui::MenuItem("Axial View", nullptr, true);
            ImGui::MenuItem("Sagittal View", nullptr, true);
            ImGui::MenuItem("Coronal View", nullptr, true);
            ImGui::EndMenu();
        }
        
        ImGui::EndMainMenuBar();
    }
}

} // namespace optirad
