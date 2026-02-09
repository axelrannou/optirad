#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <memory>

namespace optirad {

class PatientData;

class View3D {
public:
    View3D();
    ~View3D();
    
    void init();
    void render();
    void cleanup();
    
    void handleMouseInput(GLFWwindow* window);
    void handleScroll(double yOffset);
    
    void setPatientData(PatientData* data);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
    
    void initCubeRendering();
    void renderCube();
    
    // Camera parameters (orbit camera)
    float m_cameraDistance = 3.0f;
    float m_yaw = 0.0f;     // Horizontal angle in radians
    float m_pitch = 0.3f;   // Vertical angle in radians
    glm::vec3 m_target = glm::vec3(0.0f);
    
    // Mouse state
    double m_lastMouseX = 0.0;
    double m_lastMouseY = 0.0;
    bool m_leftMousePressed = false;
    bool m_rightMousePressed = false;
    bool m_middleMousePressed = false;
    double m_scrollOffset = 0.0;
    
    // Viewport size
    int m_viewportWidth = 800;
    int m_viewportHeight = 600;

    // Face labels rendering
    GLuint m_labelShaderProgram = 0;
    GLuint m_labelVAO = 0;
    GLuint m_labelVBO = 0;
    GLuint m_labelTexture = 0;
    
    void initLabelRendering();
    void renderFaceLabels();
    void createLabelTexture();
};

} // namespace optirad
