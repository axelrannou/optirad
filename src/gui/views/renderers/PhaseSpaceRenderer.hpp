#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include <mutex>

namespace optirad {

class PhaseSpaceBeamSource;

/// Visualization mode for particle coloring
enum class PhaseSpaceColorMode {
    Energy,        ///< Color by energy (blue→green→red gradient)
    ParticleType   ///< Color by type (photon=yellow, electron=blue, positron=red)
};

/// Renders phase-space particles from multiple PhaseSpaceBeamSources in the 3D viewport.
/// Supports per-beam visibility toggles (like BeamRenderer for STF beams).
/// - Particles rendered as GL_POINTS with energy-based color gradient
/// - Optional direction vectors as GL_LINES
/// - Per-beam isocenter sphere (yellow) and source sphere (orange)
class PhaseSpaceRenderer {
public:
    PhaseSpaceRenderer();
    ~PhaseSpaceRenderer();

    void init();
    void render(const glm::mat4& view, const glm::mat4& projection,
                const glm::vec3& cameraPos);
    void cleanup();

    /// Set the phase-space beam sources for rendering (one per beam).
    void setSources(const std::vector<const PhaseSpaceBeamSource*>& sources);

    // Per-beam visibility
    size_t getBeamCount() const { return m_beamVisible.size(); }
    bool isBeamVisible(size_t idx) const { return idx < m_beamVisible.size() && m_beamVisible[idx]; }
    void setBeamVisible(size_t idx, bool v) { if (idx < m_beamVisible.size()) m_beamVisible[idx] = v; m_needsRebuild = true; }
    void setAllBeamsVisible(bool v) { for (size_t i = 0; i < m_beamVisible.size(); ++i) m_beamVisible[i] = v; m_needsRebuild = true; }

    // Visualization toggles
    bool isShowingParticles() const { return m_showParticles; }
    void setShowParticles(bool v) { m_showParticles = v; }

    bool isShowingDirections() const { return m_showDirections; }
    void setShowDirections(bool v) { m_showDirections = v; m_needsRebuild = true; }

    bool isShowingAperture() const { return m_showAperture; }
    void setShowAperture(bool v) { m_showAperture = v; }

    PhaseSpaceColorMode getColorMode() const { return m_colorMode; }
    void setColorMode(PhaseSpaceColorMode mode) { m_colorMode = mode; m_needsRebuild = true; }

    float getPointSize() const { return m_pointSize; }
    void setPointSize(float s) { m_pointSize = s; }

    float getDirectionScale() const { return m_directionScale; }
    void setDirectionScale(float s) { m_directionScale = s; m_needsRebuild = true; }

    size_t getParticleCount() const { return m_totalParticleCount; }

private:
    void rebuildGeometry();
    void cleanupBeamData();

    /// Map energy to RGB color (blue→green→red gradient)
    static glm::vec3 energyToColor(float energy, float minE, float maxE);

    /// Map particle type to RGB color
    static glm::vec3 typeToColor(int particleType);

    // Shaders
    GLuint m_pointShader = 0;
    GLuint m_lineShader = 0;

    // Per-beam render data
    struct BeamRenderData {
        GLuint pointVAO = 0, pointVBO = 0;
        size_t particleCount = 0;

        GLuint dirVAO = 0, dirVBO = 0;
        size_t dirVertexCount = 0;

        glm::vec3 isocenter{0.0f};
        glm::vec3 sourcePos{0.0f};

        float gantryAngle = 0.0f;
    };
    std::vector<BeamRenderData> m_beamData;
    std::vector<bool> m_beamVisible;
    size_t m_totalParticleCount = 0;

    // Isocenter sphere mesh (shared across beams)
    GLuint m_isoVAO = 0, m_isoVBO = 0, m_isoEBO = 0;
    size_t m_isoIndexCount = 0;

    // Source pointers
    std::vector<const PhaseSpaceBeamSource*> m_sources;
    bool m_needsRebuild = false;

    // Visualization settings
    bool m_showParticles = true;
    bool m_showDirections = false;
    bool m_showAperture = false;
    PhaseSpaceColorMode m_colorMode = PhaseSpaceColorMode::Energy;
    float m_pointSize = 2.0f;
    float m_directionScale = 10.0f; // mm length of direction arrows

    mutable std::mutex m_mutex;
};

} // namespace optirad
