#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace optirad {

class AxisLabels {
public:
    AxisLabels();
    ~AxisLabels();
    
    void init();
    void render(const glm::mat4& view, const glm::mat4& projection);
    void cleanup();
    
private:
    void createLabelTexture();
    
    GLuint m_shaderProgram = 0;
    GLuint m_vao = 0;
    GLuint m_vbo = 0;
    GLuint m_texture = 0;
};

}
