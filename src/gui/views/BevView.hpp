#pragma once

#include "IView.hpp"
#include "../AppState.hpp"
#include "core/FluenceMap.hpp"
#include <GL/glew.h>
#include <imgui.h>
#include <string>
#include <vector>

namespace optirad {

/// Beam's Eye View (BEV) — displays fluence heatmap, MLC leaf positions,
/// and projected patient structure contours for a selected beam and segment.
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
    void renderBevContent();
    void updateFluenceTexture();
    void renderMlcLeaves(ImDrawList* dl, ImVec2 imgMin, ImVec2 imgMax);
    void renderProjectedContours(ImDrawList* dl, ImVec2 imgMin, ImVec2 imgMax);

    /// Jet colormap: value in [0,1] → R,G,B,A
    static void jetColormap(float t, unsigned char& r, unsigned char& g,
                            unsigned char& b);

    /// Convert BEV (X,Z) to screen coordinates.
    ImVec2 bevToScreen(double bevX, double bevZ, ImVec2 imgMin, ImVec2 imgMax) const;

    GuiAppState& m_state;

    // Selection state
    int m_beamIndex = 0;
    int m_segmentIndex = -1; // -1 = "All" (show combined fluence)

    // Zoom and pan (UV-based, same as SliceView)
    float m_zoom = 1.0f;
    float m_panU = 0.5f;
    float m_panV = 0.5f;

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

    // Contour visibility toggle
    bool m_showContours = true;

    // Leaf pair boundary cache
    std::vector<double> m_leafBoundaries; // cumulative Z positions
    int m_cachedNumLeaves = -1;

    // Dose version tracking for auto-refresh
    int m_doseVersion = -1;
};

} // namespace optirad
