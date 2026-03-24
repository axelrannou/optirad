#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "DoseColors.hpp"
#include "dose/PlanAnalysis.hpp"
#include <string>
#include <vector>

namespace optirad {

/// DVH panel: displays DVH curves for the selected dose map.
/// Supports comparison overlay when a second dose is selected.
class DVHPanel : public IPanel {
public:
    DVHPanel(GuiAppState& state);

    std::string getName() const override { return "DVH"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    // Primary dose DVH
    std::vector<StructureDoseStats> m_stats;
    std::vector<DVHCurveData> m_dvhCurves;

    // Comparison dose DVH
    std::vector<StructureDoseStats> m_compareStats;
    std::vector<DVHCurveData> m_compareDvhCurves;

    // State tracking
    int m_doseVersion = -1;
    float m_dvhMaxDose = 80.0f;

    // Per-structure visibility (int to avoid vector<bool> proxy)
    std::vector<int> m_curveVisible;

    void recompute();
    void renderDVH();
};

} // namespace optirad
