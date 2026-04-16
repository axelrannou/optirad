#include "DoseStatsPanel.hpp"
#include "../Theme.hpp"
#include <imgui.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace optirad {

DoseStatsPanel::DoseStatsPanel(GuiAppState& state) : m_state(state) {}

void DoseStatsPanel::computeStats() {
    m_stats.clear();
    m_compareStats.clear();

    auto& dm = m_state.doseManager;
    int selIdx = dm.getSelectedIdx();
    auto* sel = dm.getSelected();
    if (!sel || !sel->dose || !sel->grid || !m_state.patientData) return;

    // Primary dose stats (cached in DoseManager)
    double rxDose = (m_state.plan && m_state.plan->getPrescribedDose() > 0)
                    ? m_state.plan->getPrescribedDose() : 60.0;
    m_stats = dm.getOrComputeStats(selIdx, *m_state.patientData, rxDose);

    // Comparison dose stats (also cached)
    int cmpIdx = dm.getCompareIdx();
    auto* cmp = dm.getCompare();
    if (cmp && cmp->dose && cmp->grid) {
        m_compareStats = dm.getOrComputeStats(cmpIdx, *m_state.patientData, rxDose);
    }
}

void DoseStatsPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Dose Statistics", &m_visible);

    if (!m_state.doseAvailable()) {
        ImGui::TextDisabled("No dose data available.");
        ImGui::End();
        return;
    }

    // Auto-refresh when dose manager version changes
    int currentVersion = m_state.doseManager.version();
    if (currentVersion != m_doseVersion) {
        computeStats();
        m_doseVersion = currentVersion;
    }

    // Comparison selector
    if (m_state.doseManager.count() >= 2) {
        renderCompareSelector();
    }

    ImGui::Spacing();

    renderStatsTable();

    ImGui::End();
}

void DoseStatsPanel::renderCompareSelector() {
    auto& dm = m_state.doseManager;
    int selIdx = dm.getSelectedIdx();
    int cmpIdx = dm.getCompareIdx();

    // Build label for current comparison selection
    const char* currentLabel = "None";
    std::string cmpLabel;
    if (cmpIdx >= 0) {
        auto* entry = dm.getEntry(cmpIdx);
        if (entry) {
            cmpLabel = entry->name;
            currentLabel = cmpLabel.c_str();
        }
    }

    ImGui::Text("Compare with:");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(200.0f);
    if (ImGui::BeginCombo("##CompareSelect", currentLabel)) {
        // "None" option
        if (ImGui::Selectable("None", cmpIdx < 0)) {
            dm.setCompare(-1);
        }
        // All doses except the selected one
        for (int i = 0; i < dm.count(); ++i) {
            if (i == selIdx) continue;
            auto* entry = dm.getEntry(i);
            if (!entry) continue;
            bool isSelected = (i == cmpIdx);
            if (ImGui::Selectable(entry->name.c_str(), isSelected)) {
                dm.setCompare(i);
            }
            if (isSelected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

void DoseStatsPanel::renderStatsTable() {
    if (m_stats.empty()) {
        ImGui::TextDisabled("No statistics available.");
        return;
    }

    const auto& kColors = DoseColors::kColors;
    constexpr int kNumColors = DoseColors::kNumColors;

    bool hasComparison = !m_compareStats.empty();

    // Build a map from structure name to comparison stats for easy lookup
    auto findCompare = [&](const std::string& name) -> const StructureDoseStats* {
        for (const auto& cs : m_compareStats) {
            if (cs.name == name) return &cs;
        }
        return nullptr;
    };

    // Helper: render a colored diff value
    // For most structures: green = negative/better (lower dose), red = positive/worse (higher dose)
    // For target structures (PTV/GTV/CTV): inverted logic (higher dose is better)
    auto renderDiff = [](double primary, double compare, const char* fmt, const std::string& structName) {
        double diff = primary - compare;
        if (std::abs(diff) < 1e-6) {
            ImGui::TextDisabled("=");
        } else {
            const auto& tc = getThemeColors();
            // Check if this is a target structure (PTV, GTV, CTV)
            bool isTarget = (structName.find("PTV") != std::string::npos ||
                           structName.find("GTV") != std::string::npos ||
                           structName.find("CTV") != std::string::npos);
            
            // For targets: positive diff is good (green), negative is bad (red)
            // For OARs: negative diff is good (green), positive is bad (red)
            ImVec4 color;
            if (isTarget) {
                color = diff > 0 ? tc.passText : tc.failText;
            } else {
                color = diff > 0 ? tc.failText : tc.passText;
            }
            
            char buf[64];
            snprintf(buf, sizeof(buf), fmt, diff);
            ImGui::TextColored(color, "%s", buf);
        }
    };

    constexpr int numCols = 14;

    if (ImGui::BeginTable("DoseStatsTable", numCols,
            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
            ImGuiTableFlags_ScrollX | ImGuiTableFlags_SizingFixedFit)) {

        ImGui::TableSetupColumn("Structure",    ImGuiTableColumnFlags_None, 120.0f);
        ImGui::TableSetupColumn("Min (Gy)",     ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("Max (Gy)",     ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("Mean (Gy)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("D2% (Gy)",     ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("D5% (Gy)",     ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("D95% (Gy)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("D98% (Gy)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("V20Gy (%)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("V40Gy (%)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("V50Gy (%)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("V60Gy (%)",    ImGuiTableColumnFlags_None, 65.0f);
        ImGui::TableSetupColumn("CI",           ImGuiTableColumnFlags_None, 45.0f);
        ImGui::TableSetupColumn("HI",           ImGuiTableColumnFlags_None, 45.0f);
        ImGui::TableHeadersRow();

        for (size_t i = 0; i < m_stats.size(); ++i) {
            const auto& s = m_stats[i];
            int colorIdx = static_cast<int>(i) % kNumColors;
            const StructureDoseStats* cs = hasComparison ? findCompare(s.name) : nullptr;

            ImGui::TableNextRow();
            int col = 0;

            // Structure name
            ImGui::TableSetColumnIndex(col++);
            ImGui::TextColored(ImVec4(kColors[colorIdx][0], kColors[colorIdx][1],
                                      kColors[colorIdx][2], 1.0f),
                              "%s", s.name.c_str());

            // Helper: render value, and if comparison exists, show diff on next line
            auto cell = [&](double val, double cmpVal, const char* valFmt, const char* diffFmt) {
                ImGui::TableSetColumnIndex(col++);
                ImGui::Text(valFmt, val);
                if (cs) {
                    renderDiff(val, cmpVal, diffFmt, s.name);
                }
            };

            cell(s.minDose,  cs ? cs->minDose  : 0, "%.2f", "%+.2f");
            cell(s.maxDose,  cs ? cs->maxDose  : 0, "%.2f", "%+.2f");
            cell(s.meanDose, cs ? cs->meanDose : 0, "%.2f", "%+.2f");
            cell(s.d2,       cs ? cs->d2       : 0, "%.2f", "%+.2f");
            cell(s.d5,       cs ? cs->d5       : 0, "%.2f", "%+.2f");
            cell(s.d95,      cs ? cs->d95      : 0, "%.2f", "%+.2f");
            cell(s.d98,      cs ? cs->d98      : 0, "%.2f", "%+.2f");
            cell(s.v20,      cs ? cs->v20      : 0, "%.1f", "%+.1f");
            cell(s.v40,      cs ? cs->v40      : 0, "%.1f", "%+.1f");
            cell(s.v50,      cs ? cs->v50      : 0, "%.1f", "%+.1f");
            cell(s.v60,      cs ? cs->v60      : 0, "%.1f", "%+.1f");

            // CI
            ImGui::TableSetColumnIndex(col++);
            if (s.conformityIndex > 0.0) {
                ImGui::Text("%.3f", s.conformityIndex);
                if (cs && cs->conformityIndex > 0.0)
                    renderDiff(s.conformityIndex, cs->conformityIndex, "%+.3f", s.name);
            } else {
                ImGui::TextDisabled("-");
            }

            // HI
            ImGui::TableSetColumnIndex(col++);
            if (s.homogeneityIndex > 0.0) {
                ImGui::Text("%.3f", s.homogeneityIndex);
                if (cs && cs->homogeneityIndex > 0.0)
                    renderDiff(s.homogeneityIndex, cs->homogeneityIndex, "%+.3f", s.name);
            } else {
                ImGui::TextDisabled("-");
            }
        }

        ImGui::EndTable();
    }
}

} // namespace optirad
