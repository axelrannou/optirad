#include "DoseStatsPanel.hpp"
#include <imgui.h>
#include <algorithm>
#include <cmath>

namespace optirad {

constexpr float DoseStatsPanel::kColors[kNumColors][3];

DoseStatsPanel::DoseStatsPanel(GuiAppState& state) : m_state(state) {}

void DoseStatsPanel::computeStats() {
    m_stats.clear();
    m_dvhCurves.clear();

    if (!m_state.doseResult || !m_state.doseGrid || !m_state.patientData) return;

    // Compute comprehensive stats via PlanAnalysis
    m_stats = PlanAnalysis::computeStats(
        *m_state.doseResult, *m_state.patientData, *m_state.doseGrid);

    // Compute DVH curves
    double maxDose = m_state.doseResult->getMax();
    m_dvhMaxDose = static_cast<float>(std::ceil(maxDose / 10.0) * 10.0);
    m_dvhCurves = PlanAnalysis::computeDVHCurves(m_stats, m_dvhMaxDose);

    // Initialize visibility flags
    m_curveVisible.assign(m_dvhCurves.size(), 1);

    m_statsComputed = true;
}

void DoseStatsPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Dose Statistics & DVH", &m_visible);

    if (!m_state.optimizationDone()) {
        ImGui::TextDisabled("Run optimization first to see dose statistics.");
        ImGui::End();
        return;
    }

    // Recompute if dose data changed
    if (m_state.doseResult) {
        if (!m_statsComputed || m_state.doseResult.get() != m_lastDosePtr) {
            computeStats();
            m_lastDosePtr = m_state.doseResult.get();
        }
    }

    if (ImGui::Button("Refresh Stats")) {
        m_statsComputed = false;
        computeStats();
    }

    ImGui::Spacing();

    // Two sections: stats table and DVH
    if (ImGui::CollapsingHeader("Dose Statistics Table", ImGuiTreeNodeFlags_DefaultOpen)) {
        renderStatsTable();
    }

    if (ImGui::CollapsingHeader("DVH Curves", ImGuiTreeNodeFlags_DefaultOpen)) {
        renderDVH();
    }

    ImGui::End();
}

void DoseStatsPanel::renderStatsTable() {
    if (m_stats.empty()) {
        ImGui::TextDisabled("No statistics available.");
        return;
    }

    // Full stats table with 14 columns
    if (ImGui::BeginTable("DoseStatsTable", 14,
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

            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::TextColored(ImVec4(kColors[colorIdx][0], kColors[colorIdx][1],
                                      kColors[colorIdx][2], 1.0f),
                              "%s", s.name.c_str());

            ImGui::TableSetColumnIndex(1);  ImGui::Text("%.2f", s.minDose);
            ImGui::TableSetColumnIndex(2);  ImGui::Text("%.2f", s.maxDose);
            ImGui::TableSetColumnIndex(3);  ImGui::Text("%.2f", s.meanDose);
            ImGui::TableSetColumnIndex(4);  ImGui::Text("%.2f", s.d2);
            ImGui::TableSetColumnIndex(5);  ImGui::Text("%.2f", s.d5);
            ImGui::TableSetColumnIndex(6);  ImGui::Text("%.2f", s.d95);
            ImGui::TableSetColumnIndex(7);  ImGui::Text("%.2f", s.d98);
            ImGui::TableSetColumnIndex(8);  ImGui::Text("%.1f", s.v20);
            ImGui::TableSetColumnIndex(9);  ImGui::Text("%.1f", s.v40);
            ImGui::TableSetColumnIndex(10); ImGui::Text("%.1f", s.v50);
            ImGui::TableSetColumnIndex(11); ImGui::Text("%.1f", s.v60);
            ImGui::TableSetColumnIndex(12);
            if (s.conformityIndex > 0.0)
                ImGui::Text("%.3f", s.conformityIndex);
            else
                ImGui::TextDisabled("-");
            ImGui::TableSetColumnIndex(13);
            if (s.homogeneityIndex > 0.0)
                ImGui::Text("%.3f", s.homogeneityIndex);
            else
                ImGui::TextDisabled("-");
        }

        ImGui::EndTable();
    }
}

void DoseStatsPanel::renderDVH() {
    if (m_dvhCurves.empty()) {
        ImGui::TextDisabled("No DVH data.");
        return;
    }

    // Structure visibility toggles
    ImGui::Text("Structures:");
    for (size_t i = 0; i < m_dvhCurves.size(); ++i) {
        int colorIdx = static_cast<int>(i) % kNumColors;
        ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol_CheckMark,
            ImVec4(kColors[colorIdx][0], kColors[colorIdx][1], kColors[colorIdx][2], 1.0f));
        if (i < m_curveVisible.size()) {
            bool visible = (m_curveVisible[i] != 0);
            ImGui::Checkbox(m_dvhCurves[i].structureName.c_str(), &visible);
            m_curveVisible[i] = visible ? 1 : 0;
        }
        ImGui::PopStyleColor();
    }

    // DVH max dose slider
    ImGui::SliderFloat("Max Dose (Gy)", &m_dvhMaxDose, 10.0f, 200.0f, "%.0f");

    // DVH Canvas using ImGui draw list
    ImVec2 canvasPos = ImGui::GetCursorScreenPos();
    ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    canvasSize.y = std::max(canvasSize.y, 300.0f);
    canvasSize.x = std::max(canvasSize.x, 400.0f);

    // Margins for axis labels
    const float marginLeft = 50.0f;
    const float marginBottom = 30.0f;
    const float marginTop = 10.0f;
    const float marginRight = 10.0f;

    ImVec2 plotOrigin(canvasPos.x + marginLeft, canvasPos.y + marginTop);
    ImVec2 plotSize(canvasSize.x - marginLeft - marginRight,
                     canvasSize.y - marginTop - marginBottom);

    ImDrawList* drawList = ImGui::GetWindowDrawList();

    // Background
    drawList->AddRectFilled(
        ImVec2(plotOrigin.x, plotOrigin.y),
        ImVec2(plotOrigin.x + plotSize.x, plotOrigin.y + plotSize.y),
        IM_COL32(30, 30, 30, 255));

    // Grid lines
    for (int i = 0; i <= 10; ++i) {
        float frac = i / 10.0f;

        // Vertical grid (dose axis)
        float x = plotOrigin.x + frac * plotSize.x;
        drawList->AddLine(
            ImVec2(x, plotOrigin.y),
            ImVec2(x, plotOrigin.y + plotSize.y),
            IM_COL32(60, 60, 60, 255));

        // Dose label
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f", frac * m_dvhMaxDose);
        drawList->AddText(ImVec2(x - 10, plotOrigin.y + plotSize.y + 5), IM_COL32(200, 200, 200, 255), buf);

        // Horizontal grid (volume axis)
        float y = plotOrigin.y + plotSize.y - frac * plotSize.y;
        drawList->AddLine(
            ImVec2(plotOrigin.x, y),
            ImVec2(plotOrigin.x + plotSize.x, y),
            IM_COL32(60, 60, 60, 255));

        // Volume label
        snprintf(buf, sizeof(buf), "%.0f%%", frac * 100.0f);
        drawList->AddText(ImVec2(plotOrigin.x - 45, y - 7), IM_COL32(200, 200, 200, 255), buf);
    }

    // Border
    drawList->AddRect(
        ImVec2(plotOrigin.x, plotOrigin.y),
        ImVec2(plotOrigin.x + plotSize.x, plotOrigin.y + plotSize.y),
        IM_COL32(150, 150, 150, 255));

    // Plot DVH curves
    for (size_t ci = 0; ci < m_dvhCurves.size(); ++ci) {
        if (ci >= m_curveVisible.size() || m_curveVisible[ci] == 0) continue;

        const auto& curve = m_dvhCurves[ci];
        int colorIdx = static_cast<int>(ci) % kNumColors;
        ImU32 lineColor = IM_COL32(
            static_cast<int>(kColors[colorIdx][0] * 255),
            static_cast<int>(kColors[colorIdx][1] * 255),
            static_cast<int>(kColors[colorIdx][2] * 255),
            255);

        for (size_t j = 1; j < curve.doses.size(); ++j) {
            float x0 = plotOrigin.x + (curve.doses[j-1] / m_dvhMaxDose) * plotSize.x;
            float y0 = plotOrigin.y + plotSize.y - (curve.volumes[j-1] / 100.0f) * plotSize.y;
            float x1 = plotOrigin.x + (curve.doses[j] / m_dvhMaxDose) * plotSize.x;
            float y1 = plotOrigin.y + plotSize.y - (curve.volumes[j] / 100.0f) * plotSize.y;

            // Clip to plot area
            if (x0 <= plotOrigin.x + plotSize.x && x1 >= plotOrigin.x) {
                drawList->AddLine(ImVec2(x0, y0), ImVec2(x1, y1), lineColor, 2.0f);
            }
        }
    }

    // Reserve space for the canvas
    ImGui::Dummy(canvasSize);
}

} // namespace optirad
