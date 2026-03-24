#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <memory>
#include <mutex>
#include <vector>

namespace optirad {

class PatientData;
class Stf;
class BeamRenderer;
class PhaseSpaceRenderer;
class PhaseSpaceBeamSource;

class View3D {
public:
    View3D();
    ~View3D();
    
    void init();
    /// Render the 3D scene into an internal FBO at the given size.
    /// Call getTextureID() afterwards to display the result in ImGui::Image().
    void render(int width, int height);
    void cleanup();
    
    void handleMouseInput(GLFWwindow* window);
    void handleScroll(double yOffset);
    
    /// Handle mouse input from ImGui hover state (for docked window)
    void handleImGuiInput();
    
    void setPatientData(PatientData* data);
    void setStf(const Stf* stf);
    void setPhaseSpaceSources(const std::vector<const PhaseSpaceBeamSource*>& sources);
    BeamRenderer* getBeamRenderer();
    PhaseSpaceRenderer* getPhaseSpaceRenderer();

    /// Switch between dark and light mode (affects logo background colour)
    void setDarkMode(bool dark) { m_darkMode = dark; }
    
    /// Get the FBO color texture ID for ImGui::Image()
    GLuint getTextureID() const { return m_fboTexture; }

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
    
    void initCubeRendering();
    void renderCube();
    
    // FBO for off-screen rendering
    void initFBO();
    void resizeFBO(int width, int height);
    GLuint m_fbo = 0;
    GLuint m_fboTexture = 0;
    GLuint m_fboDepthRBO = 0;
    int m_fboWidth = 0;
    int m_fboHeight = 0;
    
    // Camera parameters (orbit camera)
    float m_cameraDistance = 3.0f;
    float m_yaw = 0.0f;     // Horizontal angle in radians
    float m_pitch = 0.3f;   // Vertical angle in radians
    glm::vec3 m_target = glm::vec3(0.0f);

    // Theme
    bool m_darkMode = true;
    
    // Mouse state - protected with mutex for thread safety
    mutable std::mutex m_mouseMutex;
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

    // Logo rendering (shown when no patient is loaded)
    void initLogoRendering();
    void renderLogo(int width, int height);
    void cleanupLogoRendering();
    GLuint m_logoShader    = 0;
    GLuint m_logoVAO       = 0;
    GLuint m_logoVBO       = 0;
    GLuint m_gradShader    = 0;
    GLuint m_gradVAO       = 0;
    GLuint m_gradVBO       = 0;
    GLuint m_logoTexDark   = 0;  // LaTIM_white_black.png
    GLuint m_logoTexLight  = 0;  // LaTIM.png
};

} // namespace optirad
