#include "StructureRenderer.hpp"
#include "core/PatientData.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>

namespace optirad {

StructureRenderer::StructureRenderer() = default;
StructureRenderer::~StructureRenderer() = default;

void StructureRenderer::init() {
    // Simple colored mesh shader
    const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        
        out vec3 FragPos;
        out vec3 Normal;
        
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        
        void main() {
            FragPos = vec3(model * vec4(aPos, 1.0));
            Normal = mat3(transpose(inverse(model))) * aNormal;
            gl_Position = projection * view * vec4(FragPos, 1.0);
        }
    )";

    const char* fragmentShaderSource = R"(
        #version 330 core
        in vec3 FragPos;
        in vec3 Normal;
        out vec4 FragColor;
        
        uniform vec3 structureColor;
        uniform vec3 lightDir;
        uniform float opacity;
        
        void main() {
            vec3 norm = normalize(Normal);
            float diff = max(dot(norm, lightDir), 0.0);
            vec3 diffuse = (0.3 + 0.7 * diff) * structureColor;
            FragColor = vec4(diffuse, opacity);
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
}

void StructureRenderer::setPatientData(PatientData* data) {
    if (m_patientData == data) return;
    m_patientData = data;
    m_needsRebuild = true;
}

void StructureRenderer::buildMeshes() {
    // Clean up old meshes
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1, &mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1, &mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1, &mesh.ebo);
    }
    m_meshes.clear();

    if (!m_patientData || !m_patientData->getStructureSet()) return;

    auto* structures = m_patientData->getStructureSet();
    for (size_t i = 0; i < structures->getCount(); ++i) {
        tessellateStructure(i);
    }

    m_needsRebuild = false;
}

void StructureRenderer::tessellateStructure(size_t structureIndex) {
    auto* structures = m_patientData->getStructureSet();
    const auto* structure = structures->getStructure(structureIndex);
    if (!structure) return;

    // Simple extrusion-based tessellation of contours
    std::vector<float> vertices;
    std::vector<unsigned int> indices;

    const auto& contours = structure->getContours();
    
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        const auto& contour = contours[ci];
        if (contour.points.size() < 3) continue;

        size_t baseIdx = vertices.size() / 6; // 6 floats per vertex (pos + normal)

        // Add vertices (convert mm to world units - divide by some scale)
        float scale = 0.001f; // 1mm = 0.001 units
        for (const auto& pt : contour.points) {
            // Position
            vertices.push_back(pt[0] * scale);
            vertices.push_back(pt[1] * scale);
            vertices.push_back(pt[2] * scale);
            
            // Normal (simple: point outward from centroid)
            vertices.push_back(0.0f);
            vertices.push_back(0.0f);
            vertices.push_back(1.0f);
        }

        // Simple triangulation (fan from first vertex)
        for (size_t i = 1; i + 1 < contour.points.size(); ++i) {
            indices.push_back(baseIdx);
            indices.push_back(baseIdx + i);
            indices.push_back(baseIdx + i + 1);
        }

        // Extrude to next contour if available (simple side walls)
        if (ci + 1 < contours.size()) {
            const auto& nextContour = contours[ci + 1];
            size_t nextBaseIdx = vertices.size() / 6;

            for (const auto& pt : nextContour.points) {
                vertices.push_back(pt[0] * scale);
                vertices.push_back(pt[1] * scale);
                vertices.push_back(pt[2] * scale);
                vertices.push_back(0.0f);
                vertices.push_back(0.0f);
                vertices.push_back(1.0f);
            }

            // Connect contours with quads
            size_t n = std::min(contour.points.size(), nextContour.points.size());
            for (size_t i = 0; i < n; ++i) {
                size_t next = (i + 1) % n;
                indices.push_back(baseIdx + i);
                indices.push_back(baseIdx + next);
                indices.push_back(nextBaseIdx + next);
                
                indices.push_back(baseIdx + i);
                indices.push_back(nextBaseIdx + next);
                indices.push_back(nextBaseIdx + i);
            }
        }
    }

    if (vertices.empty() || indices.empty()) return;

    // Create mesh
    StructureMesh mesh;
    auto color = structure->getColor();
    mesh.color = glm::vec3(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
    mesh.visible = structure->isVisible();
    mesh.indexCount = indices.size();

    glGenVertexArrays(1, &mesh.vao);
    glGenBuffers(1, &mesh.vbo);
    glGenBuffers(1, &mesh.ebo);

    glBindVertexArray(mesh.vao);
    
    glBindBuffer(GL_ARRAY_BUFFER, mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    m_meshes.push_back(mesh);
}

void StructureRenderer::render(const glm::mat4& view, const glm::mat4& projection) {
    if (!m_patientData || !m_patientData->getStructureSet()) return;

    if (m_needsRebuild) buildMeshes();
    if (m_meshes.empty()) return;

    glm::mat4 model = glm::mat4(1.0f);
    glm::vec3 lightDir = glm::normalize(glm::vec3(0.5f, 1.0f, 0.5f));

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE); // Don't write to depth for transparent objects

    glUseProgram(m_shaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniform3fv(glGetUniformLocation(m_shaderProgram, "lightDir"), 1, glm::value_ptr(lightDir));
    glUniform1f(glGetUniformLocation(m_shaderProgram, "opacity"), 0.5f);

    for (const auto& mesh : m_meshes) {
        if (!mesh.visible) continue;

        glUniform3fv(glGetUniformLocation(m_shaderProgram, "structureColor"), 1, glm::value_ptr(mesh.color));
        
        glBindVertexArray(mesh.vao);
        glDrawElements(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
}

void StructureRenderer::cleanup() {
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1, &mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1, &mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1, &mesh.ebo);
    }
    m_meshes.clear();
    
    if (m_shaderProgram) {
        glDeleteProgram(m_shaderProgram);
        m_shaderProgram = 0;
    }
}

}
