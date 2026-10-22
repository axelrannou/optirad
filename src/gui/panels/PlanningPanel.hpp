#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/workflow/DoseCalculationPipeline.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include <atomic>
#include <thread>
#include <string>
#include <vector>

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

    // Delivery technique: 0=SMLC, 1=VMAT
    int m_deliveryIdx = 0;
    struct ArcUi {
        float start = static_cast<float>(ArcSpec{}.startAngle);
        float stop = static_cast<float>(ArcSpec{}.stopAngle);
        bool clockwise = ArcSpec{}.clockwise;
        float couch = static_cast<float>(ArcSpec{}.couchAngle);
    };
    std::vector<ArcUi> m_arcs{ArcUi{}};
    float m_optSpacing = static_cast<float>(VmatSettings{}.optimizationSpacingDeg);   // degrees per fluence-optimization beam
    float m_cpSpacing = static_cast<float>(VmatSettings{}.controlPointSpacingDeg);    // degrees per control point

    // Fractions
    int m_numFractions = PlanConfig{}.numFractions;

    // Gantry angles: 0=Range, 1=List
    int m_gantryMode = 0;
    float m_gantryStart = 0.0f;
    float m_gantryStep = static_cast<float>(PlanConfig{}.gantryStep);
    float m_gantryStop = 360.0f;
    std::string m_gantryListBuf = "0 90 180 270";

    // Couch angles: 0=Range, 1=List
    int m_couchMode = 0;
    float m_couchStart = 0.0f;
    float m_couchStep = 0.0f;
    float m_couchStop = 0.0f;
    std::string m_couchListBuf = "0 0 0 0";

    // Bixel width
    float m_bixelWidth = static_cast<float>(PlanConfig{}.bixelWidth);

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
    float m_doseResolution[3] = {static_cast<float>(DoseCalcPipelineOptions{}.resolution[0]),
                                 static_cast<float>(DoseCalcPipelineOptions{}.resolution[1]),
                                 static_cast<float>(DoseCalcPipelineOptions{}.resolution[2])};
    int   m_engineIdx = 1;              // 0 = stored Dij, 1 = matrix-free
    float m_fluencePixelMm = static_cast<float>(DoseCalcPipelineOptions{}.fluencePixelMm);
    int m_computeBackendIdx = 1;
    int m_gpuMemoryBudgetMiB = 0;
    float m_relativeThreshold = static_cast<float>(100.0 * DoseCalcPipelineOptions{}.relativeThreshold);   // %
    bool  m_excludeExternal = DoseCalcPipelineOptions{}.excludeExternal;
    float m_absoluteThreshold = static_cast<float>(DoseCalcPipelineOptions{}.absoluteThreshold);
    int m_numThreads = 0;
};

} // namespace optirad
