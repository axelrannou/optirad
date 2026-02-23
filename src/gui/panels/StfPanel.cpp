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

    // ── Generate STF button ──
    if (m_isGenerating) {
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Generating STF...");
        ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));
        
        // Check if background thread finished
        if (m_generationDone) {
            m_generationThread.join();
            m_isGenerating = false;
            m_generationDone = false;
            Logger::info("STF generation complete");
        }
    } else {
        if (ImGui::Button("Generate STF", ImVec2(-1, 30))) {
            // Reset STF
            m_state.resetStf();
            m_isGenerating = true;
            m_generationDone = false;

            // Run generation on background thread to keep GUI responsive
            m_generationThread = std::thread([this]() {
                auto start = std::chrono::steady_clock::now();

                const auto& stfP = m_state.plan->getStfProperties();
                std::string radiationMode = m_state.plan->getRadiationMode();

                std::array<double, 3> iso = {0.0, 0.0, 0.0};
                if (!stfP.isoCenters.empty()) {
                    iso = stfP.isoCenters[0];
                }

                double gantryStart = !stfP.gantryAngles.empty() ? stfP.gantryAngles[0] : 0.0;
                double gantryStop = !stfP.gantryAngles.empty() ? stfP.gantryAngles.back() + 1.0 : 360.0;
                double gantryStep = stfP.gantryAngles.size() > 1
                    ? stfP.gantryAngles[1] - stfP.gantryAngles[0] : 60.0;
                double bixelWidth = stfP.bixelWidth;

                PhotonIMRTStfGenerator generator(gantryStart, gantryStep, gantryStop, bixelWidth, iso);
                generator.setMachine(m_state.plan->getMachine());
                generator.setRadiationMode(radiationMode);

                auto patientData = m_state.plan->getPatientData();
                if (patientData && patientData->hasValidCT() && patientData->hasStructures()) {
                    const auto* ct = patientData->getCTVolume();
                    const auto* structureSet = patientData->getStructureSet();
                    const auto& grid = ct->getGrid();

                    generator.setGrid(grid);
                    generator.setStructureSet(*structureSet);
                    generator.setCTResolution(grid.getSpacing());
                }

                m_state.stfProps = generator.generate();
                m_state.stf = std::make_shared<Stf>(generator.generateStf());

                auto end = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(end - start).count();
                m_statusMessage = "Generated in " + std::to_string(elapsed).substr(0, 5) + "s";

                m_generationDone = true;
            });
        }
    }

    // ── Show STF results ──
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
            if (ImGui::BeginTable("BeamTable", 4,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_ScrollY, ImVec2(0, 200))) {

                ImGui::TableSetupColumn("Beam");
                ImGui::TableSetupColumn("Gantry");
                ImGui::TableSetupColumn("Rays");
                ImGui::TableSetupColumn("Energy");
                ImGui::TableHeadersRow();

                for (size_t i = 0; i < m_state.stf->getCount(); ++i) {
                    const auto* beam = m_state.stf->getBeam(i);
                    if (!beam) continue;

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn(); ImGui::Text("%zu", i);
                    ImGui::TableNextColumn(); ImGui::Text("%.1f", beam->getGantryAngle());
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

                // Per-beam checkboxes in a scrollable region
                ImGui::BeginChild("BeamVisScroll", ImVec2(0, 200), true);
                for (size_t i = 0; i < m_beamRenderer->getBeamCount(); ++i) {
                    const auto* beam = m_state.stf->getBeam(i);
                    if (!beam) continue;

                    bool visible = m_beamRenderer->isBeamVisible(i);
                    char label[64];
                    snprintf(label, sizeof(label), "Beam %zu (%.0f deg, %zu rays)",
                             i, beam->getGantryAngle(), beam->getNumOfRays());
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
