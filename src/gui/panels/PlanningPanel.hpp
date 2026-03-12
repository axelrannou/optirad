#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include <atomic>
#include <thread>
#include <string>

namespace optirad {

/// Treatment planning panel: creates plans, generates STF, computes Dij.
///
/// Provides controls for:
///  - Radiation mode and machine selection
///  - Gantry angles (range or explicit list)
///  - Couch angles (range or explicit list)
///  - Bixel width
///  - Dose grid resolution and calculation options
class PlanningPanel : public IPanel {
public:
    PlanningPanel(GuiAppState& state);
    ~PlanningPanel();

    std::string getName() const override { return "Planning"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    // Radiation mode / machine
    int m_radiationModeIdx = 0;
    int m_machineIdx = 0;

    // Fractions
    int m_numFractions = 30;

    // Gantry angles: 0=Range, 1=List
    int m_gantryMode = 0;
    float m_gantryStart = 0.0f;
    float m_gantryStep = 4.0f;
    float m_gantryStop = 360.0f;
    std::string m_gantryListBuf = "0 90 180 270";

    // Couch angles: 0=Range, 1=List
    int m_couchMode = 0;
    float m_couchStart = 0.0f;
    float m_couchStep = 0.0f;
    float m_couchStop = 0.0f;
    std::string m_couchListBuf = "0 0 0 0";

    // Bixel width
    float m_bixelWidth = 7.0f;

    // Plan state
    bool m_planCreated = false;

    // STF generation
    bool m_isGeneratingStf = false;
    std::atomic<bool> m_stfGenerationDone{false};
    std::thread m_stfThread;
    std::string m_stfStatusMessage;

    // Dij calculation
    bool m_isCalculatingDij = false;
    std::atomic<bool> m_dijCalcDone{false};
    std::thread m_dijThread;
    std::string m_dijStatusMessage;
    int m_dijCurrentBeam = 0;
    int m_dijTotalBeams = 0;

    // Dose calc options
    float m_doseResolution[3] = {2.5f, 2.5f, 2.5f};
    float m_relativeThreshold = 1.0f;
    float m_absoluteThreshold = 1e-6f;
    int m_numThreads = 0;
};

} // namespace optirad
