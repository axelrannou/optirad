#include "DVHPanel.hpp"
#include "../Theme.hpp"
#include <imgui.h>
#include <algorithm>
#include <cmath>

namespace optirad {

DVHPanel::DVHPanel(GuiAppState& state) : m_state(state) {}

void DVHPanel::recompute() {
    m_stats.clear();
    m_dvhCurves.clear();
    m_compareStats.clear();
    m_compareDvhCurves.clear();

    auto& dm = m_state.doseStore;
    int selIdx = dm.getSelectedIdx();
    auto* sel = dm.getSelected();
    if (!sel || !sel->dose || !sel->grid || !m_state.patientData) return;

    // Primary dose (cached in GuiAppState)
    double rxDose = (m_state.plan && m_state.plan->getPrescribedDose() > 0)
                    ? m_state.plan->getPrescribedDose() : 60.0;
    m_stats = m_state.getOrComputeStats(selIdx, *m_state.patientData, rxDose);
    double maxDose = sel->dose->getMax();
    m_dvhMaxDose = static_cast<float>(maxDose);
    m_dvhCurves = PlanAnalysis::computeDVHCurves(m_stats, m_dvhMaxDose);

    // Comparison dose (also cached)
    int cmpIdx = dm.getCompareIdx();
    auto* cmp = dm.getCompare();
    if (cmp && cmp->dose && cmp->grid) {
        m_compareStats = m_state.getOrComputeStats(cmpIdx, *m_state.patientData, rxDose);
        double cmpMax = cmp->dose->getMax();
        if (cmpMax > maxDose) {
            m_dvhMaxDose = static_cast<float>(cmpMax);
        }
        m_compareDvhCurves = PlanAnalysis::computeDVHCurves(m_compareStats, m_dvhMaxDose);
    }

    // Initialize visibility flags for the max number of structures
    size_t maxStructures = std::max(m_dvhCurves.size(), m_compareDvhCurves.size());
    if (m_curveVisible.size() < maxStructures) {
        m_curveVisible.resize(maxStructures, 1);
    }
}

void DVHPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("DVH", &m_visible);

    if (!m_state.doseAvailable()) {
        ImGui::TextDisabled("No dose data available.");
        ImGui::End();
        return;
    }

    // Auto-refresh when dose manager version changes
    int currentVersion = m_state.doseStore.version();
    if (currentVersion != m_doseVersion) {
        recompute();
        m_doseVersion = currentVersion;
    }

    renderDVH();

    ImGui::End();
}

void DVHPanel::renderDVH() {
    if (m_dvhCurves.empty()) {
        ImGui::TextDisabled("No DVH data.");
        return;
    }

    const auto& kColors = DoseColors::kColors;
    constexpr int kNumColors = DoseColors::kNumColors;

    bool hasComparison = !m_compareDvhCurves.empty();

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

    if (hasComparison) {
        auto* cmpEntry = m_state.doseStore.getCompare();
        auto* selEntry = m_state.doseStore.getSelected();
        if (selEntry && cmpEntry) {
            ImGui::TextColored(ImVec4(0.6f, 0.8f, 1.0f, 1.0f),
                "Solid = %s | Dashed = %s", selEntry->name.c_str(), cmpEntry->name.c_str());
        }
    }

    // DVH max dose slider
    ImGui::SliderFloat("Max Dose (Gy)", &m_dvhMaxDose, 10.0f, 200.0f, "%.0f");

    // DVH Canvas using ImGui draw list (fills available space)
    ImVec2 canvasPos = ImGui::GetCursorScreenPos();
    ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    canvasSize.y = std::max(canvasSize.y, 200.0f);
    canvasSize.x = std::max(canvasSize.x, 300.0f);

    // Margins for axis labels
    const float marginLeft = 50.0f;
    const float marginBottom = 30.0f;
    const float marginTop = 10.0f;
    const float marginRight = 10.0f;

    ImVec2 plotOrigin(canvasPos.x + marginLeft, canvasPos.y + marginTop);
    ImVec2 plotSize(canvasSize.x - marginLeft - marginRight,
                     canvasSize.y - marginTop - marginBottom);

    ImDrawList* drawList = ImGui::GetWindowDrawList();

    const auto& tc = getThemeColors();

    // Background
    auto vecToCol32 = [](const ImVec4& v) {
        return IM_COL32(
            static_cast<int>(v.x * 255),
            static_cast<int>(v.y * 255),
            static_cast<int>(v.z * 255),
            static_cast<int>(v.w * 255));
    };

    drawList->AddRectFilled(
        ImVec2(plotOrigin.x, plotOrigin.y),
        ImVec2(plotOrigin.x + plotSize.x, plotOrigin.y + plotSize.y),
        vecToCol32(tc.dvhBackground));

    // Grid lines
    for (int i = 0; i <= 10; ++i) {
        float frac = i / 10.0f;

        // Vertical grid (dose axis)
        float x = plotOrigin.x + frac * plotSize.x;
        drawList->AddLine(
            ImVec2(x, plotOrigin.y),
            ImVec2(x, plotOrigin.y + plotSize.y),
            vecToCol32(tc.dvhGrid));

        // Dose label
        char buf[16];
        snprintf(buf, sizeof(buf), "%.0f", frac * m_dvhMaxDose);
        drawList->AddText(ImVec2(x - 10, plotOrigin.y + plotSize.y + 5), vecToCol32(tc.dvhLabel), buf);

        // Horizontal grid (volume axis)
        float y = plotOrigin.y + plotSize.y - frac * plotSize.y;
        drawList->AddLine(
            ImVec2(plotOrigin.x, y),
            ImVec2(plotOrigin.x + plotSize.x, y),
            vecToCol32(tc.dvhGrid));

        // Volume label
        snprintf(buf, sizeof(buf), "%.0f%%", frac * 100.0f);
        drawList->AddText(ImVec2(plotOrigin.x - 45, y - 7), vecToCol32(tc.dvhLabel), buf);
    }

    // Border
    drawList->AddRect(
        ImVec2(plotOrigin.x, plotOrigin.y),
        ImVec2(plotOrigin.x + plotSize.x, plotOrigin.y + plotSize.y),
        vecToCol32(tc.dvhBorder));

    // Helper lambda to draw solid DVH curves
    auto drawCurvesSolid = [&](const std::vector<DVHCurveData>& curves, float thickness) {
        for (size_t ci = 0; ci < curves.size(); ++ci) {
            if (ci >= m_curveVisible.size() || m_curveVisible[ci] == 0) continue;

            const auto& curve = curves[ci];
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

                if (x0 <= plotOrigin.x + plotSize.x && x1 >= plotOrigin.x) {
                    drawList->AddLine(ImVec2(x0, y0), ImVec2(x1, y1), lineColor, thickness);
                }
            }
        }
    };

    // Helper lambda to draw dashed DVH curves (dash=8px, gap=6px)
    auto drawCurvesDashed = [&](const std::vector<DVHCurveData>& curves, float thickness) {
        constexpr float dashLen = 8.0f;
        constexpr float gapLen  = 6.0f;
        constexpr float patternLen = dashLen + gapLen;

        for (size_t ci = 0; ci < curves.size(); ++ci) {
            if (ci >= m_curveVisible.size() || m_curveVisible[ci] == 0) continue;

            const auto& curve = curves[ci];
            int colorIdx = static_cast<int>(ci) % kNumColors;
            ImU32 lineColor = IM_COL32(
                static_cast<int>(kColors[colorIdx][0] * 255),
                static_cast<int>(kColors[colorIdx][1] * 255),
                static_cast<int>(kColors[colorIdx][2] * 255),
                255);

            float accumLen = 0.0f; // running length along the polyline for dash pattern

            for (size_t j = 1; j < curve.doses.size(); ++j) {
                float x0 = plotOrigin.x + (curve.doses[j-1] / m_dvhMaxDose) * plotSize.x;
                float y0 = plotOrigin.y + plotSize.y - (curve.volumes[j-1] / 100.0f) * plotSize.y;
                float x1 = plotOrigin.x + (curve.doses[j] / m_dvhMaxDose) * plotSize.x;
                float y1 = plotOrigin.y + plotSize.y - (curve.volumes[j] / 100.0f) * plotSize.y;

                float dx = x1 - x0, dy = y1 - y0;
                float segLen = std::sqrt(dx * dx + dy * dy);
                if (segLen < 0.5f) { accumLen += segLen; continue; }
                float ux = dx / segLen, uy = dy / segLen;

                float t = 0.0f;
                while (t < segLen) {
                    float posInPattern = std::fmod(accumLen + t, patternLen);
                    bool inDash = (posInPattern < dashLen);
                    float remaining = inDash ? (dashLen - posInPattern) : (patternLen - posInPattern);
                    float tEnd = std::min(t + remaining, segLen);

                    if (inDash) {
                        float sx = x0 + ux * t,  sy = y0 + uy * t;
                        float ex = x0 + ux * tEnd, ey = y0 + uy * tEnd;
                        if (sx <= plotOrigin.x + plotSize.x && ex >= plotOrigin.x) {
                            drawList->AddLine(ImVec2(sx, sy), ImVec2(ex, ey), lineColor, thickness);
                        }
                    }
                    t = tEnd;
                }
                accumLen += segLen;
            }
        }
    };

    // Draw comparison curves first (dashed, behind)
    if (hasComparison) {
        drawCurvesDashed(m_compareDvhCurves, 1.5f);
    }

    // Draw primary curves (solid, thicker)
    drawCurvesSolid(m_dvhCurves, 2.0f);

    // Reserve space for the canvas
    ImGui::Dummy(canvasSize);
}

} // namespace optirad
