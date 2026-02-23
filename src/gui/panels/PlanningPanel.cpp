#include "PlanningPanel.hpp"
#include "io/MachineLoader.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>
#include <chrono>

namespace optirad {

static const char* kRadiationModes[] = { "photons" };
static const char* kMachines[] = { "Generic", "Varian_TrueBeam6MV" };
static const bool kIsPhaseSpace[] = { false, true };

PlanningPanel::PlanningPanel(GuiAppState& state) : m_state(state) {}

void PlanningPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Planning", &m_visible);

    // Gate: need DICOM data first
    if (!m_state.dicomLoaded()) {
        ImGui::TextDisabled("Load DICOM data first to configure a plan.");
        ImGui::End();
        return;
    }

    // ── Radiation Mode ──
    if (ImGui::CollapsingHeader("Radiation Mode", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        ImGui::Combo("Mode", &m_radiationModeIdx, kRadiationModes, IM_ARRAYSIZE(kRadiationModes));
        ImGui::Combo("Machine", &m_machineIdx, kMachines, IM_ARRAYSIZE(kMachines));
        if (kIsPhaseSpace[m_machineIdx]) {
            ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f), "Phase-space (IAEA PSF)");
        } else {
            ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.5f, 1.0f), "Generic (pencil-beam)");
        }
        ImGui::Unindent();
    }

    // ── Fractions ──
    if (ImGui::CollapsingHeader("Fractions", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        ImGui::InputInt("Num Fractions", &m_numFractions);
        if (m_numFractions < 1) m_numFractions = 1;
        ImGui::Unindent();
    }

    // ── Gantry Angles ──
    if (ImGui::CollapsingHeader("Gantry Angles", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        ImGui::DragFloat("Start (deg)", &m_gantryStart, 1.0f, 0.0f, 360.0f, "%.1f");
        ImGui::DragFloat("Step (deg)", &m_gantryStep, 0.5f, 0.5f, 90.0f, "%.1f");
        ImGui::DragFloat("Stop (deg)", &m_gantryStop, 1.0f, 1.0f, 720.0f, "%.1f");

        // Show resulting number of beams
        int numBeams = 0;
        if (m_gantryStep > 0.0f && m_gantryStop > m_gantryStart) {
            numBeams = static_cast<int>(std::ceil((m_gantryStop - m_gantryStart) / m_gantryStep));
        }
        ImGui::Text("Number of beams: %d", numBeams);
        ImGui::Unindent();
    }

    // ── Bixel Width (only for generic machines) ──
    if (!kIsPhaseSpace[m_machineIdx]) {
        if (ImGui::CollapsingHeader("Bixel Width", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Indent();
            ImGui::DragFloat("Width (mm)", &m_bixelWidth, 0.5f, 1.0f, 20.0f, "%.1f");
            ImGui::Unindent();
        }
    }

    ImGui::Separator();
    ImGui::Spacing();

    // ── Validation ──
    bool paramsValid = (m_gantryStep > 0.0f) &&
                       (m_gantryStop > m_gantryStart) &&
                       (kIsPhaseSpace[m_machineIdx] || m_bixelWidth > 0.0f) &&
                       (m_numFractions > 0);

    if (!paramsValid) {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "Fix parameter errors above.");
    }

    // ── Create Plan button ──
    if (!paramsValid) ImGui::BeginDisabled();

    if (ImGui::Button("Create Plan", ImVec2(-1, 30))) {
        // Reset downstream state
        m_state.resetPlan();

        auto plan = std::make_shared<Plan>();
        plan->setName("TreatmentPlan");
        plan->setRadiationMode(kRadiationModes[m_radiationModeIdx]);
        plan->setNumOfFractions(m_numFractions);
        plan->setPatientData(m_state.patientData);

        // Load machine from JSON data file
        try {
            plan->setMachine(MachineLoader::load(
                kRadiationModes[m_radiationModeIdx],
                kMachines[m_machineIdx]));
        } catch (const std::exception& e) {
            Logger::error("Failed to load machine: " + std::string(e.what()));
            ImGui::End();
            return;
        }

        // STF properties
        StfProperties stfProps;
        stfProps.setGantryAngles(m_gantryStart, m_gantryStep, m_gantryStop);
        stfProps.bixelWidth = m_bixelWidth;

        // Compute isocenter from target structures
        auto iso = plan->computeIsoCenter();
        stfProps.setUniformIsoCenter(iso);

        plan->setStfProperties(stfProps);

        m_state.plan = plan;
        m_planCreated = true;

        Logger::info("Plan created: " + std::to_string(static_cast<int>(
            std::ceil((m_gantryStop - m_gantryStart) / m_gantryStep))) +
            " beams, bixelWidth=" + std::to_string(m_bixelWidth) + "mm");
    }

    if (!paramsValid) ImGui::EndDisabled();

    // Show plan status
    if (m_state.planCreated()) {
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "Plan created.");

        const auto& stfProps = m_state.plan->getStfProperties();
        ImGui::Text("  Beams: %zu", stfProps.gantryAngles.size());
        ImGui::Text("  Iso: (%.1f, %.1f, %.1f) mm",
                    stfProps.isoCenters[0][0],
                    stfProps.isoCenters[0][1],
                    stfProps.isoCenters[0][2]);
    }

    // ── Generate STF ──
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("STF Generation", ImGuiTreeNodeFlags_DefaultOpen)) {
        bool canGenerateStf = m_state.planCreated() &&
                              !m_state.isPhaseSpaceMachine();

        if (!m_state.planCreated()) {
            ImGui::TextDisabled("Create a plan first to generate STF.");
        } else if (m_state.isPhaseSpaceMachine()) {
            ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f),
                "Phase-space machine: use the Phase Space panel instead.");
        }

        if (m_isGeneratingStf) {
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Generating STF...");
            ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));

            if (m_stfGenerationDone) {
                m_stfThread.join();
                m_isGeneratingStf = false;
                m_stfGenerationDone = false;
                Logger::info("STF generation complete");
            }
        } else {
            if (!canGenerateStf) ImGui::BeginDisabled();

            if (ImGui::Button("Generate STF", ImVec2(-1, 30))) {
                m_state.resetStf();
                m_isGeneratingStf = true;
                m_stfGenerationDone = false;

                m_stfThread = std::thread([this]() {
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
                    m_stfStatusMessage = "Generated in " + std::to_string(elapsed).substr(0, 5) + "s";

                    m_stfGenerationDone = true;
                });
            }

            if (!canGenerateStf) ImGui::EndDisabled();
        }

        // Show STF status
        if (m_state.stfGenerated()) {
            if (!m_stfStatusMessage.empty()) {
                ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "%s", m_stfStatusMessage.c_str());
            }
            ImGui::Text("Beams: %zu | Rays: %zu | Bixels: %zu",
                        m_state.stf->getCount(),
                        m_state.stf->getTotalNumOfRays(),
                        m_state.stf->getTotalNumOfBixels());
        }
    }

    ImGui::End();
}

} // namespace optirad
