#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "DoseColors.hpp"
#include "core/workflow/PlanAnalysis.hpp"
#include <string>
#include <vector>

namespace optirad {

/// Dose statistics panel: shows comprehensive plan analysis after optimization.
/// Supports comparison of two dose maps side-by-side.
class DoseStatsPanel : public IPanel {
public:
    DoseStatsPanel(GuiAppState& state);

    std::string getName() const override { return "Dose Statistics"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    // Primary dose stats
    std::vector<StructureDoseStats> m_stats;

    // Comparison dose stats
    std::vector<StructureDoseStats> m_compareStats;

    // State tracking for auto-refresh
    int m_doseVersion = -1;

    void computeStats();
    void renderStatsTable();
    void renderCompareSelector();
};

} // namespace optirad
