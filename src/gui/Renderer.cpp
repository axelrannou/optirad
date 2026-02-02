#include "Renderer.hpp"
#include "Logger.hpp"

namespace optirad {

Renderer::Renderer() = default;
Renderer::~Renderer() { shutdown(); }

bool Renderer::init() {
    // TODO: Initialize OpenGL context
    // TODO: Initialize Dear ImGui
    // ImGui::CreateContext();
    // ImGui_ImplGlfw_InitForOpenGL(window, true);
    // ImGui_ImplOpenGL3_Init("#version 330");
    
    m_initialized = true;
    Logger::info("Renderer initialized");
    return true;
}

void Renderer::shutdown() {
    if (m_initialized) {
        // TODO: ImGui_ImplOpenGL3_Shutdown();
        // TODO: ImGui_ImplGlfw_Shutdown();
        // TODO: ImGui::DestroyContext();
        m_initialized = false;
    }
}

void Renderer::beginFrame() {
    // TODO: glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // TODO: ImGui_ImplOpenGL3_NewFrame();
    // TODO: ImGui_ImplGlfw_NewFrame();
    // TODO: ImGui::NewFrame();
}

void Renderer::endFrame() {
    // TODO: ImGui::Render();
    // TODO: ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

} // namespace optirad
