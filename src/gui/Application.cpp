#include "Application.hpp"
#include "panels/PatientPanel.hpp"
#include "panels/BeamPanel.hpp"
#include "panels/OptimizationPanel.hpp"
#include "panels/DoseStatsPanel.hpp"
#include "panels/LogPanel.hpp"
#include "Logger.hpp"

namespace optirad {

Application::Application() = default;
Application::~Application() = default;

bool Application::init() {
    Logger::info("Initializing OptiRad application...");
    
    m_window = std::make_unique<Window>();
    if (!m_window->create("OptiRad", 1920, 1080)) {
        Logger::error("Failed to create window");
        return false;
    }
    
    m_renderer = std::make_unique<Renderer>();
    if (!m_renderer->init()) {
        Logger::error("Failed to initialize renderer");
        return false;
    }
    
    // Create panels
    m_panels.push_back(std::make_unique<PatientPanel>());
    m_panels.push_back(std::make_unique<BeamPanel>());
    m_panels.push_back(std::make_unique<OptimizationPanel>());
    m_panels.push_back(std::make_unique<DoseStatsPanel>());
    m_panels.push_back(std::make_unique<LogPanel>());
    
    m_running = true;
    Logger::info("Application initialized successfully");
    return true;
}

void Application::run() {
    while (m_running && !m_window->shouldClose()) {
        m_window->pollEvents();
        update();
        render();
        m_window->swapBuffers();
    }
}

void Application::update() {
    for (auto& panel : m_panels) {
        panel->update();
    }
}

void Application::render() {
    m_renderer->beginFrame();
    
    for (auto& panel : m_panels) {
        panel->render();
    }
    
    m_renderer->endFrame();
}

void Application::shutdown() {
    Logger::info("Shutting down application...");
    m_panels.clear();
    m_renderer.reset();
    m_window.reset();
}

} // namespace optirad
