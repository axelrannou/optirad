#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace optirad {

class CubeRenderer {
public:
    CubeRenderer();
    ~CubeRenderer();
    
    void init();
    void render(const glm::mat4& view, const glm::mat4& projection);
    void cleanup();
    
private:
    GLuint m_shaderProgram = 0;
    GLuint m_vao = 0;
    GLuint m_vbo = 0;
    GLuint m_ebo = 0;
};

}
