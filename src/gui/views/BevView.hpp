#pragma once

#include "IView.hpp"
#include "../AppState.hpp"
#include "core/FluenceMap.hpp"
#include "renderers/BevRenderer.hpp"
#include <GL/glew.h>
#include <imgui.h>
#include <string>
#include <vector>

namespace optirad {

/// Beam's Eye View with source-perspective anatomy/MLC and a fluence display.
class BevView : public IView {
public:
    explicit BevView(GuiAppState& state);
    ~BevView() override;

    void render() override;
    void update() override {}
    void resize(int, int) override {}
    std::string getName() const override { return "BEV"; }

private:
    void renderControls();
    void renderArcMotion();
    void renderSceneControls();
    void handleZoom(float& zoom);
    void renderScene(const Beam& beam, const std::vector<const Aperture*>& apertures,
                     const std::vector<double>& boundaries, double extentX, double extentZ);
    bool m_arcPlaying = false;
    float m_arcTime = 0.0f;
    float m_arcPlaybackSpeed = 1.0f;
    int m_arcSelection = 0;
    int m_selectedDoseId = -1;
    void renderBevContent();
    void updateFluenceTexture();
    void renderMlcLeaves(ImDrawList* dl, ImVec2 imgMin, ImVec2 imgMax);

    /// Jet colormap: value in [0,1] → R,G,B,A
    static void jetColormap(float t, unsigned char& r, unsigned char& g,
                            unsigned char& b);

    /// Convert BEV (X,Z) to screen coordinates.
    ImVec2 bevToScreen(double bevX, double bevZ, ImVec2 imgMin, ImVec2 imgMax) const;

    GuiAppState& m_state;

    // Selection state
    int m_beamIndex = 0;
    int m_segmentIndex = 0; // -1 = "All" (show combined fluence/aperture outlines)

    // Centered fluence zoom
    float m_zoom = 1.0f;

    // BEV extent for the current beam (mm)
    double m_bevMinX = 0, m_bevMaxX = 0;
    double m_bevMinZ = 0, m_bevMaxZ = 0;

    // Fluence texture
    GLuint m_fluenceTexID = 0;
    int m_texWidth = 0;
    int m_texHeight = 0;
    int m_lastTexWidth = -1;
    int m_lastTexHeight = -1;
    bool m_needsUpdate = true;

    BevRenderer m_sceneRenderer;
    bool m_showAnatomy = true;
    float m_sceneZoom = 1.0f;

    // Dose version tracking for auto-refresh
    int m_doseVersion = -1;
};

} // namespace optirad
