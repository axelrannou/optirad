#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>

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
    void buildMeshes();
    void tessellateStructure(size_t structureIndex);
    
    struct StructureMesh {
        GLuint vao = 0;
        GLuint vbo = 0;
        GLuint ebo = 0;
        size_t indexCount = 0;
        glm::vec3 color;
        bool visible = true;
    };
    
    PatientData* m_patientData = nullptr;
    bool m_needsRebuild = false;
    
    std::vector<StructureMesh> m_meshes;
    GLuint m_shaderProgram = 0;
};

}
