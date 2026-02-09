#include "View3D.hpp"
#include "renderers/CubeRenderer.hpp"
#include "renderers/AxisLabels.hpp"
#include <imgui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <cmath>
#include <memory>

namespace optirad {

// Use pimpl pattern to hide renderer details
struct View3D::Impl {
    std::unique_ptr<CubeRenderer> cubeRenderer;
    std::unique_ptr<AxisLabels> axisLabels;
};

View3D::View3D() : m_impl(std::make_unique<Impl>()) {}
View3D::~View3D() = default;

void View3D::init() {
    glEnable(GL_DEPTH_TEST);
    glLineWidth(2.0f);
    
    m_impl->cubeRenderer = std::make_unique<CubeRenderer>();
    m_impl->cubeRenderer->init();
    
    m_impl->axisLabels = std::make_unique<AxisLabels>();
    m_impl->axisLabels->init();
}

void View3D::handleScroll(double yOffset) {
    m_scrollOffset = yOffset;
}

void View3D::handleMouseInput(GLFWwindow* window) {
    glfwGetFramebufferSize(window, &m_viewportWidth, &m_viewportHeight);
    
    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);
    
    bool leftPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    bool rightPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    bool middlePressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
    
    if (ImGui::GetIO().WantCaptureMouse) {
        m_leftMousePressed = false;
        m_rightMousePressed = false;
        m_middleMousePressed = false;
        m_scrollOffset = 0.0;
        return;
    }
    
    double deltaX = mouseX - m_lastMouseX;
    double deltaY = mouseY - m_lastMouseY;
    
    if (leftPressed && m_leftMousePressed) {
        const float sensitivity = 0.005f;
        m_yaw -= static_cast<float>(deltaX) * sensitivity;
        m_pitch += static_cast<float>(deltaY) * sensitivity;
        const float maxPitch = glm::half_pi<float>() - 0.01f;
        m_pitch = std::clamp(m_pitch, -maxPitch, maxPitch);
    }
    
    if (rightPressed && m_rightMousePressed) {
        const float panSpeed = 0.002f * m_cameraDistance;
        glm::vec3 forward(std::sin(m_yaw), 0.0f, std::cos(m_yaw));
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        m_target -= right * static_cast<float>(deltaX) * panSpeed;
        m_target += up * static_cast<float>(deltaY) * panSpeed;
    }
    
    if (m_scrollOffset != 0.0) {
        const float zoomSpeed = 0.1f;
        m_cameraDistance -= static_cast<float>(m_scrollOffset) * zoomSpeed * m_cameraDistance;
        m_cameraDistance = std::clamp(m_cameraDistance, 0.5f, 50.0f);
        m_scrollOffset = 0.0;
    }
    
    m_lastMouseX = mouseX;
    m_lastMouseY = mouseY;
    m_leftMousePressed = leftPressed;
    m_rightMousePressed = rightPressed;
    m_middleMousePressed = middlePressed;
}

void View3D::render() {
    // Calculate camera matrices
    float cosPitch = std::cos(m_pitch), sinPitch = std::sin(m_pitch);
    float cosYaw = std::cos(m_yaw), sinYaw = std::sin(m_yaw);
    
    glm::vec3 cameraOffset(
        m_cameraDistance * cosPitch * sinYaw,
        m_cameraDistance * sinPitch,
        m_cameraDistance * cosPitch * cosYaw
    );
    
    glm::vec3 cameraPos = m_target + cameraOffset;
    glm::mat4 view = glm::lookAt(cameraPos, m_target, glm::vec3(0.0f, 1.0f, 0.0f));
    
    float aspect = static_cast<float>(m_viewportWidth) / static_cast<float>(std::max(m_viewportHeight, 1));
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);
    
    // Render components
    m_impl->cubeRenderer->render(view, projection);
    m_impl->axisLabels->render(view, projection);
}

void View3D::cleanup() {
    m_impl->cubeRenderer->cleanup();
    m_impl->axisLabels->cleanup();
}

}
