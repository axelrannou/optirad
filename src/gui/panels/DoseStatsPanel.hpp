#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "dose/PlanAnalysis.hpp"
#include <string>
#include <vector>

namespace optirad {

/// Dose statistics panel: shows comprehensive plan analysis and DVH curves after optimization.
class DoseStatsPanel : public IPanel {
public:
    DoseStatsPanel(GuiAppState& state);

    std::string getName() const override { return "Dose Statistics"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    std::vector<StructureDoseStats> m_stats;
    std::vector<DVHCurveData> m_dvhCurves;
    bool m_statsComputed = false;
    const DoseMatrix* m_lastDosePtr = nullptr;  // to detect when dose data changes

    // DVH display settings
    float m_dvhMaxDose = 80.0f;
    bool m_showDVH = true;

    // Per-structure visibility for DVH (int instead of bool to avoid vector<bool> proxy)
    std::vector<int> m_curveVisible;

    void computeStats();
    void renderStatsTable();
    void renderDVH();

    // Structure colors (cycle through these)
    static constexpr int kNumColors = 8;
    static constexpr float kColors[kNumColors][3] = {
        {1.0f, 0.0f, 0.0f},  // red
        {0.0f, 0.7f, 0.0f},  // green
        {0.2f, 0.4f, 1.0f},  // blue
        {1.0f, 0.8f, 0.0f},  // yellow
        {1.0f, 0.4f, 0.0f},  // orange
        {0.8f, 0.0f, 0.8f},  // purple
        {0.0f, 0.8f, 0.8f},  // cyan
        {0.6f, 0.3f, 0.1f},  // brown
    };
};

} // namespace optirad
