#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/Plan.hpp"
#include "core/Machine.hpp"
#include <array>
#include <string>
#include <atomic>
#include <thread>

namespace optirad {

/// Planning panel: configure treatment plan parameters after DICOM is loaded.
/// Shows radiation mode, machine, gantry angles, bixel width, fractions.
/// "Create Plan" button activates once all parameters are valid.
/// "Generate STF" button activates once the plan is created (non-phase-space).
class PlanningPanel : public IPanel {
public:
    PlanningPanel(GuiAppState& state);

    std::string getName() const override { return "Planning"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    // Plan parameters (editable)
    int m_radiationModeIdx = 0;  // 0 = photons
    int m_machineIdx = 0;        // 0 = Generic
    int m_numFractions = 30;
    float m_gantryStart = 0.0f;
    float m_gantryStep = 4.0f;
    float m_gantryStop = 360.0f;
    float m_bixelWidth = 7.0f;

    // Derived
    bool m_planCreated = false;

    // STF generation state
    std::atomic<bool> m_isGeneratingStf{false};
    std::atomic<bool> m_stfGenerationDone{false};
    std::thread m_stfThread;
    std::string m_stfStatusMessage;
};

} // namespace optirad
