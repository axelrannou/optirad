#include "Application.hpp"
#include "panels/PatientPanel.hpp"
#include "views/SliceView.hpp"
#include "utils/Logger.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <GL/gl.h>

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
    m_window = glfwCreateWindow(1920, 1080, "OptiRad TPS", nullptr, nullptr);
    if (!m_window) {
        Logger::error("Failed to create GLFW window");
        glfwTerminate();
        return false;
    }
    
    glfwMakeContextCurrent(m_window);
    glfwSwapInterval(1); // VSync
    
    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    
    ImGui_ImplGlfw_InitForOpenGL(m_window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    
    // Create panels
    m_patientPanel = std::make_unique<PatientPanel>();
    
    // Create slice views
    m_axialView = std::make_unique<SliceView>(SliceOrientation::Axial);
    m_sagittalView = std::make_unique<SliceView>(SliceOrientation::Sagittal);
    m_coronalView = std::make_unique<SliceView>(SliceOrientation::Coronal);
    
    return true;
}

void Application::run() {
    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();
        
        // Start ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // Render main menu
        renderMenuBar();
        
        // Render panels
        m_patientPanel->render();
        
        // Update views if patient data changed
        if (m_patientPanel->hasData()) {
            auto* data = m_patientPanel->getPatientData();
            m_axialView->setPatientData(data);
            m_sagittalView->setPatientData(data);
            m_coronalView->setPatientData(data);
        }
        
        // Render views
        m_axialView->render();
        m_sagittalView->render();
        m_coronalView->render();
        
        // Render ImGui
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(m_window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        glfwSwapBuffers(m_window);
    }
}

void Application::shutdown() {
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
