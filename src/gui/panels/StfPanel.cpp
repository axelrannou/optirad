#include "StfPanel.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "../views/renderers/BeamRenderer.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <chrono>

namespace optirad {

StfPanel::StfPanel(GuiAppState& state) : m_state(state) {}

void StfPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("STF Generation", &m_visible);

    // Gate: need plan first
    if (!m_state.planCreated()) {
        ImGui::TextDisabled("Create a plan first in the Planning panel.");
        ImGui::End();
        return;
    }

    // Gate: STF generation only for generic machines
    if (m_state.isPhaseSpaceMachine()) {
        ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f),
            "Machine uses phase-space beam source.");
        ImGui::TextWrapped(
            "STF generation is not needed for phase-space machines. "
            "Use the Phase Space panel to load and inspect beam data.");
        ImGui::End();
        return;
    }

    // Show plan summary
    const auto& plan = m_state.plan;
    const auto& stfProps = plan->getStfProperties();
    ImGui::Text("Plan: %s", plan->getName().c_str());
    ImGui::Text("Mode: %s | Beams: %zu | Bixel: %.1f mm",
                plan->getRadiationMode().c_str(),
                stfProps.gantryAngles.size(),
                stfProps.bixelWidth);

    ImGui::Separator();
    ImGui::Spacing();

    // ── Show STF results ──
    if (!m_state.stfGenerated()) {
        ImGui::TextDisabled("Generate STF from the Planning panel.");
    }
    if (m_state.stfGenerated()) {
        ImGui::Spacing();
        ImGui::Separator();
        
        if (!m_statusMessage.empty()) {
            ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "%s", m_statusMessage.c_str());
        }

        ImGui::Text("Beams:  %zu", m_state.stf->getCount());
        ImGui::Text("Rays:   %zu", m_state.stf->getTotalNumOfRays());
        ImGui::Text("Bixels: %zu", m_state.stf->getTotalNumOfBixels());

        ImGui::Spacing();

        // Per-beam details in a scrollable table
        if (ImGui::CollapsingHeader("Beam Details")) {
            if (ImGui::BeginTable("BeamTable", 5,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_ScrollY, ImVec2(0, 0))) {

                ImGui::TableSetupColumn("Beam");
                ImGui::TableSetupColumn("Gantry");
                ImGui::TableSetupColumn("Couch");
                ImGui::TableSetupColumn("Rays");
                ImGui::TableSetupColumn("Energy");
                ImGui::TableHeadersRow();

                for (size_t i = 0; i < m_state.stf->getCount(); ++i) {
                    const auto* beam = m_state.stf->getBeam(i);
                    if (!beam) continue;

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::Text("%zu", i);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", beam->getGantryAngle());
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", beam->getCouchAngle());
                    ImGui::TableNextColumn(); ImGui::Text("%zu", beam->getNumOfRays());
                    ImGui::TableNextColumn();
                    if (beam->getNumOfRays() > 0) {
                        ImGui::Text("%.1f MV", beam->getRay(0)->getEnergy());
                    }
                }
                ImGui::EndTable();
            }
        }


        // ── Beam visibility controls ──
        if (m_beamRenderer && m_beamRenderer->getBeamCount() > 0) {
            ImGui::Spacing();
            if (ImGui::CollapsingHeader("Beam Visibility", ImGuiTreeNodeFlags_DefaultOpen)) {
                // Select All / Deselect All buttons
                if (ImGui::Button("Show All")) {
                    m_beamRenderer->setAllBeamsVisible(true);
                }
                ImGui::SameLine();
                if (ImGui::Button("Hide All")) {
                    m_beamRenderer->setAllBeamsVisible(false);
                }

                ImGui::Separator();

                // Per-beam checkboxes in a scrollable region (fill remaining space)
                float remainH = ImGui::GetContentRegionAvail().y;
                ImGui::BeginChild("BeamVisScroll", ImVec2(-1, remainH > 60.0f ? remainH : 200.0f), true);
                for (size_t i = 0; i < m_beamRenderer->getBeamCount(); ++i) {
                    const auto* beam = m_state.stf->getBeam(i);
                    if (!beam) continue;

                    bool visible = m_beamRenderer->isBeamVisible(i);
                    char label[80];
                    snprintf(label, sizeof(label), "Beam %zu (G%.0f C%.0f, %zu rays)",
                             i, beam->getGantryAngle(), beam->getCouchAngle(), beam->getNumOfRays());
                    if (ImGui::Checkbox(label, &visible)) {
                        m_beamRenderer->setBeamVisible(i, visible);
                    }
                }
                ImGui::EndChild();
            }
        }
    }

    ImGui::End();
}

} // namespace optirad
