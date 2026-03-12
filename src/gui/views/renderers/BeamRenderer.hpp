#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include <mutex>

namespace optirad {

class Stf;

/// Renders beams and rays from an Stf in the 3D view.
/// - Isocenter: yellow sphere
/// - Per-beam rays: green lines from source → each ray position
/// - Per-beam visibility toggles
class BeamRenderer {
public:
    BeamRenderer();
    ~BeamRenderer();

    void init();
    void render(const glm::mat4& view, const glm::mat4& projection,
                const glm::vec3& cameraPos);
    void cleanup();

    void setStf(const Stf* stf);

    // Per-beam visibility
    size_t getBeamCount() const { return m_beamVisible.size(); }
    bool isBeamVisible(size_t idx) const { return idx < m_beamVisible.size() && m_beamVisible[idx]; }
    void setBeamVisible(size_t idx, bool v) { if (idx < m_beamVisible.size()) m_beamVisible[idx] = v; }
    void setAllBeamsVisible(bool v) { for (size_t i = 0; i < m_beamVisible.size(); ++i) m_beamVisible[i] = v; }

private:
    void rebuildGeometry();

    // Shaders
    GLuint m_lineShader = 0;

    // Isocenter sphere
    GLuint m_isoVAO = 0, m_isoVBO = 0, m_isoEBO = 0;
    size_t m_isoIndexCount = 0;
    glm::vec3 m_isocenter{0.0f};

    // Per-beam ray data (source → each ray position)
    struct BeamRayData {
        GLuint vao = 0, vbo = 0;
        size_t vertexCount = 0;
        glm::vec3 sourcePos{0.0f};
        float gantryAngle = 0.0f;
        float couchAngle = 0.0f;
        size_t numRays = 0;
    };
    std::vector<BeamRayData> m_beamRays;
    std::vector<bool> m_beamVisible;

    const Stf* m_stf = nullptr;
    bool m_needsRebuild = false;

    mutable std::mutex m_mutex;
};

} // namespace optirad
