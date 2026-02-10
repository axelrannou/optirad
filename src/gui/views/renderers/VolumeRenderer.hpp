#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <cstdint>

namespace optirad {

class PatientData;

class VolumeRenderer {
public:
    VolumeRenderer();
    ~VolumeRenderer();

    void init();
    void render(const glm::mat4& view, const glm::mat4& projection, const glm::vec3& cameraPos);
    void cleanup();

    void setPatientData(PatientData* data);

    // Window/level for transfer function
    void setWindowLevel(int windowWidth, int windowCenter);

private:
    void uploadVolumeTexture();
    void createProxyCube();
    void createShaders();

    PatientData* m_patientData = nullptr;
    bool m_needsUpload = false;

    // 3D texture
    GLuint m_volumeTexture = 0;

    // Proxy cube geometry (front and back faces for raycasting)
    GLuint m_cubeVAO = 0;
    GLuint m_cubeVBO = 0;
    GLuint m_cubeEBO = 0;

    // Framebuffer for back-face positions
    GLuint m_backFBO = 0;
    GLuint m_backTexture = 0;
    GLuint m_backDepthRBO = 0;

    // Shaders
    GLuint m_positionShader = 0;  // Renders cube positions to texture
    GLuint m_raycastShader = 0;   // Raycasts through volume

    // Render quad for final pass
    GLuint m_quadVAO = 0;
    GLuint m_quadVBO = 0;

    // Volume dimensions
    int m_dimX = 0, m_dimY = 0, m_dimZ = 0;
    glm::vec3 m_volumeScale = glm::vec3(1.0f);

    // Window/level - better defaults for CT visualization
    int m_windowWidth = 1500;  // Changed from 400
    int m_windowCenter = -400; // Changed from 40 (good for lung tissue)

    // Viewport for FBO resize
    int m_lastWidth = 0, m_lastHeight = 0;

    void createBackFBO(int width, int height);
};

}
