#include "Window.hpp"
#include "Logger.hpp"

namespace optirad {

Window::Window() = default;
Window::~Window() { destroy(); }

bool Window::create(const std::string& title, int width, int height) {
    m_width = width;
    m_height = height;
    
    // TODO: Initialize GLFW and create window
    // glfwInit();
    // m_handle = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    
    Logger::info("Window created: " + title + " (" + 
                 std::to_string(width) + "x" + std::to_string(height) + ")");
    return true;
}

void Window::destroy() {
    // TODO: glfwDestroyWindow and glfwTerminate
    m_handle = nullptr;
}

void Window::pollEvents() {
    // TODO: glfwPollEvents();
}

void Window::swapBuffers() {
    // TODO: glfwSwapBuffers(m_handle);
}

bool Window::shouldClose() const {
    // TODO: return glfwWindowShouldClose(m_handle);
    return false;
}

int Window::getWidth() const { return m_width; }
int Window::getHeight() const { return m_height; }

} // namespace optirad
