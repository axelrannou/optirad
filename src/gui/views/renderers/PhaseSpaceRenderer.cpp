#include "PhaseSpaceRenderer.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "phsp/PhaseSpaceData.hpp"
#include "utils/Logger.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace optirad {

// mm → GL world units (same as BeamRenderer / StructureRenderer)
static constexpr float kScale = 0.001f;

PhaseSpaceRenderer::PhaseSpaceRenderer() = default;
PhaseSpaceRenderer::~PhaseSpaceRenderer() = default;

// ─── Sphere mesh generation (same as BeamRenderer) ──────────────────
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
            verts.push_back(x); verts.push_back(y); verts.push_back(z);
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

// ─── Shader compilation ─────────────────────────────────────────────
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

// ─── Color mapping ──────────────────────────────────────────────────
glm::vec3 PhaseSpaceRenderer::energyToColor(float energy, float minE, float maxE) {
    if (maxE <= minE) return glm::vec3(0.5f, 0.5f, 0.5f);
    float t = std::clamp((energy - minE) / (maxE - minE), 0.0f, 1.0f);

    // Blue(0) → Cyan(0.25) → Green(0.5) → Yellow(0.75) → Red(1.0)
    glm::vec3 color;
    if (t < 0.25f) {
        float s = t / 0.25f;
        color = glm::vec3(0.0f, s, 1.0f);
    } else if (t < 0.5f) {
        float s = (t - 0.25f) / 0.25f;
        color = glm::vec3(0.0f, 1.0f, 1.0f - s);
    } else if (t < 0.75f) {
        float s = (t - 0.5f) / 0.25f;
        color = glm::vec3(s, 1.0f, 0.0f);
    } else {
        float s = (t - 0.75f) / 0.25f;
        color = glm::vec3(1.0f, 1.0f - s, 0.0f);
    }
    return color;
}

glm::vec3 PhaseSpaceRenderer::typeToColor(int particleType) {
    switch (particleType) {
        case 1:  return glm::vec3(1.0f, 1.0f, 0.3f);
        case -1: return glm::vec3(0.3f, 0.5f, 1.0f);
        case 2:  return glm::vec3(1.0f, 0.3f, 0.3f);
        default: return glm::vec3(0.7f, 0.7f, 0.7f);
    }
}

// ─── Init ───────────────────────────────────────────────────────────
void PhaseSpaceRenderer::init() {
    const char* pointVert = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aColor;
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform float pointSize;
        out vec3 vColor;
        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
            gl_PointSize = pointSize;
            vColor = aColor;
        }
    )";
    const char* pointFrag = R"(
        #version 330 core
        in vec3 vColor;
        out vec4 FragColor;
        void main() {
            vec2 coord = gl_PointCoord - vec2(0.5);
            if (dot(coord, coord) > 0.25) discard;
            FragColor = vec4(vColor, 0.8);
        }
    )";

    m_pointShader = compileShader(pointVert, pointFrag);

    const char* lineVert = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aColor;
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        out vec3 vColor;
        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
            vColor = aColor;
        }
    )";
    const char* lineFrag = R"(
        #version 330 core
        in vec3 vColor;
        uniform float alpha;
        out vec4 FragColor;
        void main() { FragColor = vec4(vColor, alpha); }
    )";

    m_lineShader = compileShader(lineVert, lineFrag);

    // Build shared isocenter sphere mesh
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

    glEnable(GL_PROGRAM_POINT_SIZE);
}

// ─── Set sources ────────────────────────────────────────────────────
void PhaseSpaceRenderer::setSources(const std::vector<const PhaseSpaceBeamSource*>& sources) {
    std::lock_guard<std::mutex> lock(m_mutex);

    // Check if sources actually changed
    if (sources.size() == m_sources.size()) {
        bool same = true;
        for (size_t i = 0; i < sources.size(); ++i) {
            if (sources[i] != m_sources[i]) { same = false; break; }
        }
        if (same) return;
    }

    m_sources = sources;
    m_beamVisible.resize(sources.size(), true);
    m_needsRebuild = true;
}

// ─── Cleanup per-beam GL data ───────────────────────────────────────
void PhaseSpaceRenderer::cleanupBeamData() {
    for (auto& bd : m_beamData) {
        if (bd.pointVAO) { glDeleteVertexArrays(1, &bd.pointVAO); bd.pointVAO = 0; }
        if (bd.pointVBO) { glDeleteBuffers(1, &bd.pointVBO); bd.pointVBO = 0; }
        if (bd.dirVAO) { glDeleteVertexArrays(1, &bd.dirVAO); bd.dirVAO = 0; }
        if (bd.dirVBO) { glDeleteBuffers(1, &bd.dirVBO); bd.dirVBO = 0; }
    }
    m_beamData.clear();
    m_totalParticleCount = 0;
}

// ─── Build geometry from phase-space data (all beams) ───────────────
void PhaseSpaceRenderer::rebuildGeometry() {
    cleanupBeamData();

    m_beamData.resize(m_sources.size());
    m_totalParticleCount = 0;

    for (size_t b = 0; b < m_sources.size(); ++b) {
        const auto* source = m_sources[b];
        if (!source || !source->isBuilt()) continue;

        auto& bd = m_beamData[b];

        const auto& vizData = source->getVisualizationSample();
        const auto& particles = vizData.particles();
        if (particles.empty()) continue;

        // Store isocenter and source position
        const auto& iso = source->getIsocenter();
        bd.isocenter = glm::vec3(
            static_cast<float>(iso[0]) * kScale,
            static_cast<float>(iso[1]) * kScale,
            static_cast<float>(iso[2]) * kScale);

        const auto& sp = source->getSourcePosition();
        bd.sourcePos = glm::vec3(
            static_cast<float>(sp[0]) * kScale,
            static_cast<float>(sp[1]) * kScale,
            static_cast<float>(sp[2]) * kScale);

        bd.gantryAngle = static_cast<float>(source->getGantryAngle());

        // Compute energy range for color mapping
        const auto& metrics = source->getMetrics();
        float minE = static_cast<float>(metrics.minEnergy);
        float maxE = static_cast<float>(metrics.maxEnergy);

        // ── Build particle points (pos3 + color3 interleaved) ──
        std::vector<float> pointData;
        pointData.reserve(particles.size() * 6);

        for (const auto& p : particles) {
            pointData.push_back(static_cast<float>(p.position[0]) * kScale);
            pointData.push_back(static_cast<float>(p.position[1]) * kScale);
            pointData.push_back(static_cast<float>(p.position[2]) * kScale);

            glm::vec3 color;
            if (m_colorMode == PhaseSpaceColorMode::Energy) {
                color = energyToColor(static_cast<float>(p.energy), minE, maxE);
            } else {
                color = typeToColor(static_cast<int>(p.type));
            }
            pointData.push_back(color.r);
            pointData.push_back(color.g);
            pointData.push_back(color.b);
        }

        glGenVertexArrays(1, &bd.pointVAO);
        glGenBuffers(1, &bd.pointVBO);
        glBindVertexArray(bd.pointVAO);
        glBindBuffer(GL_ARRAY_BUFFER, bd.pointVBO);
        glBufferData(GL_ARRAY_BUFFER, pointData.size() * sizeof(float),
                     pointData.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
                              reinterpret_cast<void*>(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glBindVertexArray(0);

        bd.particleCount = particles.size();
        m_totalParticleCount += particles.size();

        // ── Build direction lines (if enabled) ──
        if (m_showDirections) {
            std::vector<float> dirData;
            dirData.reserve(particles.size() * 12);

            float scale = m_directionScale * kScale;

            for (const auto& p : particles) {
                float px = static_cast<float>(p.position[0]) * kScale;
                float py = static_cast<float>(p.position[1]) * kScale;
                float pz = static_cast<float>(p.position[2]) * kScale;

                float dx = static_cast<float>(p.direction[0]) * scale;
                float dy = static_cast<float>(p.direction[1]) * scale;
                float dz = static_cast<float>(p.direction[2]) * scale;

                glm::vec3 color;
                if (m_colorMode == PhaseSpaceColorMode::Energy) {
                    color = energyToColor(static_cast<float>(p.energy), minE, maxE);
                } else {
                    color = typeToColor(static_cast<int>(p.type));
                }

                dirData.push_back(px); dirData.push_back(py); dirData.push_back(pz);
                dirData.push_back(color.r); dirData.push_back(color.g); dirData.push_back(color.b);

                dirData.push_back(px + dx); dirData.push_back(py + dy); dirData.push_back(pz + dz);
                dirData.push_back(color.r * 0.5f); dirData.push_back(color.g * 0.5f); dirData.push_back(color.b * 0.5f);
            }

            glGenVertexArrays(1, &bd.dirVAO);
            glGenBuffers(1, &bd.dirVBO);
            glBindVertexArray(bd.dirVAO);
            glBindBuffer(GL_ARRAY_BUFFER, bd.dirVBO);
            glBufferData(GL_ARRAY_BUFFER, dirData.size() * sizeof(float),
                         dirData.data(), GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
                                  reinterpret_cast<void*>(3 * sizeof(float)));
            glEnableVertexAttribArray(1);
            glBindVertexArray(0);

            bd.dirVertexCount = particles.size() * 2;
        }
    }

    m_needsRebuild = false;
    Logger::info("PhaseSpaceRenderer: built geometry for " +
                 std::to_string(m_sources.size()) + " beams, " +
                 std::to_string(m_totalParticleCount) + " particles total");
}

// ─── Render ─────────────────────────────────────────────────────────
void PhaseSpaceRenderer::render(const glm::mat4& view, const glm::mat4& projection,
                                 const glm::vec3& cameraPos) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_sources.empty()) return;
    if (m_needsRebuild) rebuildGeometry();
    if (m_beamData.empty()) return;

    // Same -90° X rotation as other renderers for LPS→OpenGL
    glm::mat4 model = glm::rotate(glm::mat4(1.0f), -glm::half_pi<float>(),
                                   glm::vec3(1.0f, 0.0f, 0.0f));

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (size_t b = 0; b < m_beamData.size(); ++b) {
        if (b >= m_beamVisible.size() || !m_beamVisible[b]) continue;

        const auto& bd = m_beamData[b];
        if (bd.particleCount == 0) continue;

        // ── Isocenter sphere (yellow) ──
        {
            glUseProgram(m_lineShader);
            float sphereRadius = 0.008f;
            glm::mat4 isoModel = model *
                glm::translate(glm::mat4(1.0f), bd.isocenter) *
                glm::scale(glm::mat4(1.0f), glm::vec3(sphereRadius));

            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                               1, GL_FALSE, glm::value_ptr(isoModel));
            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "view"),
                               1, GL_FALSE, glm::value_ptr(view));
            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "projection"),
                               1, GL_FALSE, glm::value_ptr(projection));
            glUniform1f(glGetUniformLocation(m_lineShader, "alpha"), 1.0f);

            glVertexAttrib3f(1, 1.0f, 1.0f, 0.0f); // yellow

            glBindVertexArray(m_isoVAO);
            glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_isoIndexCount),
                           GL_UNSIGNED_INT, nullptr);
        }

        // ── Source point sphere (orange) ──
        {
            float sphereRadius = 0.006f;
            glm::mat4 srcModel = model *
                glm::translate(glm::mat4(1.0f), bd.sourcePos) *
                glm::scale(glm::mat4(1.0f), glm::vec3(sphereRadius));

            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                               1, GL_FALSE, glm::value_ptr(srcModel));
            glVertexAttrib3f(1, 1.0f, 0.5f, 0.0f); // orange

            glBindVertexArray(m_isoVAO);
            glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_isoIndexCount),
                           GL_UNSIGNED_INT, nullptr);
        }

        // ── Particle points ──
        if (m_showParticles && bd.pointVAO) {
            glUseProgram(m_pointShader);
            glUniformMatrix4fv(glGetUniformLocation(m_pointShader, "model"),
                               1, GL_FALSE, glm::value_ptr(model));
            glUniformMatrix4fv(glGetUniformLocation(m_pointShader, "view"),
                               1, GL_FALSE, glm::value_ptr(view));
            glUniformMatrix4fv(glGetUniformLocation(m_pointShader, "projection"),
                               1, GL_FALSE, glm::value_ptr(projection));
            glUniform1f(glGetUniformLocation(m_pointShader, "pointSize"), m_pointSize);

            glBindVertexArray(bd.pointVAO);
            glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(bd.particleCount));
        }

        // ── Direction vectors ──
        if (m_showDirections && bd.dirVAO && bd.dirVertexCount > 0) {
            glUseProgram(m_lineShader);
            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "model"),
                               1, GL_FALSE, glm::value_ptr(model));
            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "view"),
                               1, GL_FALSE, glm::value_ptr(view));
            glUniformMatrix4fv(glGetUniformLocation(m_lineShader, "projection"),
                               1, GL_FALSE, glm::value_ptr(projection));
            glUniform1f(glGetUniformLocation(m_lineShader, "alpha"), 0.4f);

            glLineWidth(1.0f);
            glBindVertexArray(bd.dirVAO);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(bd.dirVertexCount));
        }
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_BLEND);
}

// ─── Cleanup ────────────────────────────────────────────────────────
void PhaseSpaceRenderer::cleanup() {
    cleanupBeamData();
    if (m_isoVAO) { glDeleteVertexArrays(1, &m_isoVAO); m_isoVAO = 0; }
    if (m_isoVBO) { glDeleteBuffers(1, &m_isoVBO); m_isoVBO = 0; }
    if (m_isoEBO) { glDeleteBuffers(1, &m_isoEBO); m_isoEBO = 0; }
    if (m_pointShader) { glDeleteProgram(m_pointShader); m_pointShader = 0; }
    if (m_lineShader) { glDeleteProgram(m_lineShader); m_lineShader = 0; }
}

} // namespace optirad
