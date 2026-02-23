#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include <atomic>
#include <thread>
#include <string>

namespace optirad {

class PhaseSpaceRenderer;

/// Phase-space beam source panel: loads and inspects IAEA phase-space data.
/// Generates one PhaseSpaceBeamSource per gantry angle defined in the plan
/// (mirroring how generateStf creates one beam per gantry angle).
///
/// Provides controls for:
///  - Collimator & couch angles (shared across all beams)
///  - Max particles & viz sample per beam
///  - Per-beam visibility toggles (show/hide like STF beams)
///  - Aggregate particle statistics and energy histogram
///  - Verification metrics
class PhaseSpacePanel : public IPanel {
public:
    PhaseSpacePanel(GuiAppState& state);

    std::string getName() const override { return "Phase Space"; }
    void update() override {}
    void render() override;

    void setPhaseSpaceRenderer(PhaseSpaceRenderer* renderer) { m_phaseSpaceRenderer = renderer; }

private:
    void renderLoadControls();
    void renderBeamVisibility();
    void renderVisualizationControls();
    void renderStatistics();
    void renderEnergyHistogram();

    GuiAppState& m_state;
    PhaseSpaceRenderer* m_phaseSpaceRenderer = nullptr;

    // Loading state
    std::atomic<bool> m_isLoading{false};
    std::atomic<bool> m_loadDone{false};
    std::thread m_loadThread;
    std::string m_statusMessage;

    // Configuration (collimator/couch shared across all beams)
    float m_collimatorAngle = 0.0f;
    float m_couchAngle = 0.0f;
    int m_maxParticlesK = 1000;    // thousands (per beam)
    int m_vizSampleSizeK = 100;    // thousands (per beam)

    // Currently selected beam for stats display (-1 = aggregate)
    int m_selectedBeam = -1;
};

} // namespace optirad
