#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include <mutex>

namespace optirad {

class PatientData;

class StructureRenderer {
public:
    StructureRenderer();
    ~StructureRenderer();
    
    void init();
    void render(const glm::mat4& view, const glm::mat4& projection);
    void cleanup();
    
    void setPatientData(PatientData* data);

private:
    void buildMeshes();              // Public wrapper with mutex lock
    void buildMeshes_unlocked();     // Internal implementation (requires m_meshMutex held)
    void tessellateStructure(size_t structureIndex);
    
    struct StructureMesh {
        GLuint vao = 0;
        GLuint vbo = 0;
        GLuint ebo = 0;
        size_t indexCount = 0;
        glm::vec3 color;
        bool visible = true;
        glm::vec3 center;
    };

    // WBOIT framebuffer
    GLuint m_wboitFBO    = 0;
    GLuint m_accumTex    = 0;
    GLuint m_revealTex   = 0;
    GLuint m_depthRBO    = 0;
    int    m_wboitWidth  = 0;
    int    m_wboitHeight = 0;

    // Shaders
    GLuint m_accumProgram     = 0;  // replaces m_shaderProgram
    GLuint m_compositeProgram = 0;

    // Full-screen quad
    GLuint m_quadVao = 0;
    GLuint m_quadVbo = 0;

    // Cached uniform locations
    GLint m_uModel, m_uView, m_uProjection, m_uViewPos;
    GLint m_uNormalMatrix, m_uLightDir1, m_uLightDir2;
    GLint m_uColor, m_uOpacity;

    // New private methods
    void ensureWboitFBO(int width, int height);
    void destroyWboitFBO();
    
    PatientData* m_patientData = nullptr;
    bool m_needsRebuild = false;
    
    // Thread safety for mesh vector accessed from multiple threads
    mutable std::mutex m_meshMutex;
    std::vector<StructureMesh> m_meshes;
    GLuint m_shaderProgram = 0;
};

}
