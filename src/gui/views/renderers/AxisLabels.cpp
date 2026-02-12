#include "AxisLabels.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <vector>

namespace optirad {

// Bitmap font data for S, I, P, A, R, L (8x8 pixels each)
static const unsigned char g_letterBitmaps[6][8] = {
    {0b01111110, 0b11000000, 0b11000000, 0b01111100, 0b00000110, 0b00000110, 0b11111100, 0b00000000}, // S
    {0b11111111, 0b00011000, 0b00011000, 0b00011000, 0b00011000, 0b00011000, 0b11111111, 0b00000000}, // I
    {0b11111100, 0b11000110, 0b11000110, 0b11111100, 0b11000000, 0b11000000, 0b11000000, 0b00000000}, // P
    {0b00111100, 0b01100110, 0b11000011, 0b11111111, 0b11000011, 0b11000011, 0b11000011, 0b00000000}, // A
    {0b11111100, 0b11000110, 0b11000110, 0b11111100, 0b11001100, 0b11000110, 0b11000011, 0b00000000}, // R
    {0b11000000, 0b11000000, 0b11000000, 0b11000000, 0b11000000, 0b11000000, 0b11111111, 0b00000000}, // L
};

AxisLabels::AxisLabels() = default;
AxisLabels::~AxisLabels() = default;

void AxisLabels::createLabelTexture() {
    const int letterWidth = 8, letterHeight = 8;
    const int atlasWidth = 6 * letterWidth, atlasHeight = letterHeight;
    
    std::vector<unsigned char> textureData(atlasWidth * atlasHeight * 4, 0);
    
    for (int letter = 0; letter < 6; ++letter) {
        for (int y = 0; y < letterHeight; ++y) {
            unsigned char row = g_letterBitmaps[letter][y];
            for (int x = 0; x < letterWidth; ++x) {
                bool pixel = (row >> (7 - x)) & 1;
                int idx = (y * atlasWidth + letter * letterWidth + x) * 4;
                if (pixel) {
                    textureData[idx + 0] = textureData[idx + 1] = textureData[idx + 2] = textureData[idx + 3] = 255;
                } else {
                    textureData[idx + 3] = 0;
                }
            }
        }
    }
    
    glGenTextures(1, &m_texture);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, atlasWidth, atlasHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void AxisLabels::init() {
    createLabelTexture();
    
    const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
            TexCoord = aTexCoord;
        }
    )";
    
    const char* fragmentShaderSource = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D labelTexture;
        uniform vec3 labelColor;
        void main() {
            vec4 texColor = texture(labelTexture, TexCoord);
            if (texColor.a < 0.1) discard;
            FragColor = vec4(labelColor, texColor.a);
        }
    )";
    
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertexShaderSource, nullptr);
    glCompileShader(vs);
    
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fs);
    
    m_shaderProgram = glCreateProgram();
    glAttachShader(m_shaderProgram, vs);
    glAttachShader(m_shaderProgram, fs);
    glLinkProgram(m_shaderProgram);
    glDeleteShader(vs);
    glDeleteShader(fs);
    
    float quadVertices[] = {
        -0.3f, -0.3f, 0.0f,  0.0f, 1.0f,
         0.3f, -0.3f, 0.0f,  1.0f, 1.0f,
         0.3f,  0.3f, 0.0f,  1.0f, 0.0f,
        -0.3f, -0.3f, 0.0f,  0.0f, 1.0f,
         0.3f,  0.3f, 0.0f,  1.0f, 0.0f,
        -0.3f,  0.3f, 0.0f,  0.0f, 0.0f,
    };
    
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
    
    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(0);
}

void AxisLabels::render(const glm::mat4& view, const glm::mat4& projection) {
    struct Label { glm::vec3 pos; int idx; glm::vec3 color; };
    const float offset = 0.60f;
    Label labels[] = {
        {{offset, 0.0f, 0.0f}, 5, {1.0f, 0.5f, 0.5f}}, // L (was R)
        {{-offset, 0.0f, 0.0f}, 4, {0.5f, 1.0f, 0.5f}}, // R (was L)
        {{0.0f, offset, 0.0f}, 0, {0.5f, 0.5f, 1.0f}}, // S
        {{0.0f, -offset, 0.0f}, 1, {1.0f, 1.0f, 0.5f}}, // I
        {{0.0f, 0.0f, offset}, 3, {1.0f, 0.5f, 1.0f}}, // A
        {{0.0f, 0.0f, -offset}, 2, {0.5f, 1.0f, 1.0f}}, // P
    };
    
    glm::vec3 camRight(view[0][0], view[1][0], view[2][0]);
    glm::vec3 camUp(view[0][1], view[1][1], view[2][1]);
    const float labelSize = 0.02f;
    
    glUseProgram(m_shaderProgram);
    glBindVertexArray(m_vao);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    
    for (const auto& label : labels) {
        glm::mat4 model(1.0f);
        model[0] = glm::vec4(camRight * labelSize, 0.0f);
        model[1] = glm::vec4(camUp * labelSize, 0.0f);
        model[2] = glm::vec4(glm::normalize(glm::cross(camRight, camUp)) * labelSize, 0.0f);
        model[3] = glm::vec4(label.pos, 1.0f);
        
        glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(m_shaderProgram, "labelColor"), 1, glm::value_ptr(label.color));
        
        float u0 = label.idx / 6.0f, u1 = (label.idx + 1) / 6.0f;
        float quadVertices[] = {
            -1.0f, -1.0f, 0.0f,  u0, 1.0f,
             1.0f, -1.0f, 0.0f,  u1, 1.0f,
             1.0f,  1.0f, 0.0f,  u1, 0.0f,
            -1.0f, -1.0f, 0.0f,  u0, 1.0f,
             1.0f,  1.0f, 0.0f,  u1, 0.0f,
            -1.0f,  1.0f, 0.0f,  u0, 0.0f,
        };
        
        // Allocate sufficient VBO size or use sub-data safely
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(quadVertices), quadVertices);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }
    
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glBindVertexArray(0);
    glUseProgram(0);
}

void AxisLabels::cleanup() {
    glDeleteVertexArrays(1, &m_vao);
    glDeleteBuffers(1, &m_vbo);
    glDeleteTextures(1, &m_texture);
    glDeleteProgram(m_shaderProgram);
}

}
