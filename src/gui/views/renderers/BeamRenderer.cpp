#include "BeamRenderer.hpp"
#include "core/Stf.hpp"
#include "core/Beam.hpp"
#include "core/Ray.hpp"
#include "utils/Logger.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace optirad {

// mm → GL world units (same as StructureRenderer)
static constexpr float kScale = 0.001f;

BeamRenderer::BeamRenderer() = default;
BeamRenderer::~BeamRenderer() = default;

// ─── Sphere mesh generation ─────────────────────────────────────────
static void generateSphere(std::vector<float>& verts, std::vector<unsigned int>& idxs,
                           int slices = 16, int stacks = 12) {
    for (int i = 0; i <= stacks; ++i) {
        float phi = static_cast<float>(M_PI) * i / stacks;
        float y = std::cos(phi);
        float r = std::sin(phi);
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * static_cast<float>(M_PI) * j / slices;
            float x = r * std::cos(theta);
            float z = r * std::sin(theta);
            verts.push_back(x);
            verts.push_back(y);
            verts.push_back(z);
        }
    }
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            int a = i * (slices + 1) + j;
            int b = a + slices + 1;
            idxs.push_back(a); idxs.push_back(b); idxs.push_back(a + 1);
            idxs.push_back(a + 1); idxs.push_back(b); idxs.push_back(b + 1);
        }
    }
}

// ─── Shader compilation helper ──────────────────────────────────────
static GLuint compileShader(const char* vertSrc, const char* fragSrc) {
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertSrc, nullptr);
    glCompileShader(vs);

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragSrc, nullptr);
    glCompileShader(fs);

    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

// ─── Init ───────────────────────────────────────────────────────────
void BeamRenderer::init() {
    // Simple line / point shader: uniform color, model/view/proj
    const char* lineVert = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
        }
    )";
    const char* lineFrag = R"(
        #version 330 core
        uniform vec4 uColor;
        out vec4 FragColor;
        void main() { FragColor = uColor; }
    )";

    m_lineShader = compileShader(lineVert, lineFrag);

    // Build a unit sphere for isocenter / source points
    std::vector<float> sphereVerts;
    std::vector<unsigned int> sphereIdxs;
    generateSphere(sphereVerts, sphereIdxs);

    glGenVertexArrays(1, &m_isoVAO);
    glGenBuffers(1, &m_isoVBO);
    glGenBuffers(1, &m_isoEBO);

    glBindVertexArray(m_isoVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_isoVBO);
    glBufferData(GL_ARRAY_BUFFER, sphereVerts.size() * sizeof(float),
                 sphereVerts.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_isoEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sphereIdxs.size() * sizeof(unsigned int),
                 sphereIdxs.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    m_isoIndexCount = sphereIdxs.size();
}

// ─── Set STF data ───────────────────────────────────────────────────
void BeamRenderer::setStf(const Stf* stf) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_stf == stf) return;
    m_stf = stf;
    m_needsRebuild = true;
}

// ─── Build geometry from STF ────────────────────────────────────────
void BeamRenderer::rebuildGeometry() {
    // Cleanup old beam data
    for (auto& bd : m_beamRays) {
        if (bd.vao) glDeleteVertexArrays(1, &bd.vao);
        if (bd.vbo) glDeleteBuffers(1, &bd.vbo);
    }
    m_beamRays.clear();
    m_beamVisible.clear();

    if (!m_stf || m_stf->isEmpty()) {
        m_needsRebuild = false;
        return;
    }

    const auto& beams = m_stf->getBeams();
    m_isocenter = glm::vec3(
        static_cast<float>(beams[0].getIsocenter()[0]) * kScale,
        static_cast<float>(beams[0].getIsocenter()[1]) * kScale,
        static_cast<float>(beams[0].getIsocenter()[2]) * kScale
    );

    // ── Per-beam: actual ray lines (source → each ray position in LPS) ──
    for (size_t bi = 0; bi < beams.size(); ++bi) {
        const auto& beam = beams[bi];
        const auto& rays = beam.getRays();
        if (rays.empty()) continue;

        auto src = beam.getSourcePoint();
        glm::vec3 sourceGL(static_cast<float>(src[0]) * kScale,
                            static_cast<float>(src[1]) * kScale,
                            static_cast<float>(src[2]) * kScale);

        // Build line segments: source → rayPos for each ray
        std::vector<float> lineVerts;
        lineVerts.reserve(rays.size() * 6); // 2 vertices × 3 floats per ray

        for (const auto& ray : rays) {
            auto rp = ray.getRayPos(); // LPS coords in mm
            glm::vec3 rpGL(static_cast<float>(rp[0]) * kScale,
                            static_cast<float>(rp[1]) * kScale,
                            static_cast<float>(rp[2]) * kScale);

            lineVerts.push_back(sourceGL.x);
            lineVerts.push_back(sourceGL.y);
            lineVerts.push_back(sourceGL.z);
            lineVerts.push_back(rpGL.x);
            lineVerts.push_back(rpGL.y);
            lineVerts.push_back(rpGL.z);
        }

        BeamRayData bd;
        bd.sourcePos = sourceGL;
        bd.gantryAngle = static_cast<float>(beam.getGantryAngle());
        bd.couchAngle = static_cast<float>(beam.getCouchAngle());
        bd.vertexCount = lineVerts.size() / 3;
        bd.numRays = rays.size();

        glGenVertexArrays(1, &bd.vao);
        glGenBuffers(1, &bd.vbo);
        glBindVertexArray(bd.vao);
        glBindBuffer(GL_ARRAY_BUFFER, bd.vbo);
        glBufferData(GL_ARRAY_BUFFER, lineVerts.size() * sizeof(float),
                     lineVerts.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);

        m_beamRays.push_back(bd);
    }

    // All beams visible by default
    m_beamVisible.assign(m_beamRays.size(), true);

    m_needsRebuild = false;
    Logger::info("BeamRenderer: built geometry for " + std::to_string(beams.size()) +
                 " beams (" + std::to_string(m_stf->getTotalNumOfRays()) + " rays)");
}

// ─── Render ─────────────────────────────────────────────────────────
void BeamRenderer::render(const glm::mat4& view, const glm::mat4& projection,
                          const glm::vec3& cameraPos) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (!m_stf || m_stf->isEmpty()) return;
    if (m_needsRebuild) rebuildGeometry();
    if (m_beamRays.empty()) return;

    glUseProgram(m_lineShader);

    // Same -90° X rotation as StructureRenderer for consistent coordinate system
    glm::mat4 model = glm::rotate(glm::mat4(1.0f), -glm::half_pi<float>(),
                                   glm::vec3(1.0f, 0.0f, 0.0f));

    glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                       1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // ── Isocenter sphere (yellow) ──
    {
        float sphereRadius = 0.008f; // world units (= 8mm)
        glm::mat4 isoModel = model *
            glm::translate(glm::mat4(1.0f), m_isocenter) *
            glm::scale(glm::mat4(1.0f), glm::vec3(sphereRadius));

        glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                           1, GL_FALSE, glm::value_ptr(isoModel));
        glUniform4f(glGetUniformLocation(m_lineShader, "uColor"),
                    1.0f, 1.0f, 0.0f, 1.0f); // yellow

        glBindVertexArray(m_isoVAO);
        glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_isoIndexCount),
                       GL_UNSIGNED_INT, nullptr);
    }

    // Reset model for lines
    glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                       1, GL_FALSE, glm::value_ptr(model));

    // ── Per-beam ray lines (green) ──
    glLineWidth(1.0f);
    size_t totalBeams = m_beamRays.size();
    for (size_t i = 0; i < totalBeams; ++i) {
        if (i < m_beamVisible.size() && !m_beamVisible[i]) continue;

        const auto& bd = m_beamRays[i];
        if (bd.vertexCount == 0) continue;

        glUniform4f(glGetUniformLocation(m_lineShader, "uColor"),
                    0.3f, 1.0f, 0.3f, 0.5f); // green, semi-transparent

        glBindVertexArray(bd.vao);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(bd.vertexCount));
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_BLEND);
    glLineWidth(1.0f);
}

// ─── Cleanup ────────────────────────────────────────────────────────
void BeamRenderer::cleanup() {
    if (m_isoVAO) { glDeleteVertexArrays(1, &m_isoVAO); m_isoVAO = 0; }
    if (m_isoVBO) { glDeleteBuffers(1, &m_isoVBO); m_isoVBO = 0; }
    if (m_isoEBO) { glDeleteBuffers(1, &m_isoEBO); m_isoEBO = 0; }

    for (auto& bd : m_beamRays) {
        if (bd.vao) glDeleteVertexArrays(1, &bd.vao);
        if (bd.vbo) glDeleteBuffers(1, &bd.vbo);
    }
    m_beamRays.clear();

    if (m_lineShader) { glDeleteProgram(m_lineShader); m_lineShader = 0; }
}

} // namespace optirad
