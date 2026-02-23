#include "PlanningPanel.hpp"
#include "io/MachineLoader.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>

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

    ImGui::End();
}

} // namespace optirad
