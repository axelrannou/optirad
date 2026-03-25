#include "OptimizationPanel.hpp"
#include "Theme.hpp"
#include "optimization/ObjectiveBuilder.hpp"
#include "geometry/Structure.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>
#include <sstream>

namespace optirad {

static const char* kObjectiveTypes[] = {
    "Squared Deviation", "Squared Overdose", "Squared Underdose", "MinDVH", "MaxDVH"
};

OptimizationPanel::OptimizationPanel(GuiAppState& state) : m_state(state) {}

OptimizationPanel::~OptimizationPanel() {
    if (m_optThread.joinable()) m_optThread.join();
}

void OptimizationPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Optimization", &m_visible);

    if (!m_state.dijComputed()) {
        ImGui::TextDisabled("Calculate Dij first to enable optimization.");
        ImGui::End();
        return;
    }

    // ── Initialize objectives from curated structure list ──
    // Only include clinically relevant structures with protocol-appropriate defaults.
    // User can still add/remove objectives via the table UI.
    if (!m_objectivesInitialized && m_state.patientData && m_state.patientData->getStructureSet()) {
        const auto* ss = m_state.patientData->getStructureSet();
        
        // Curated default objectives: structure name pattern -> {typeIdx, dose, weight, volPct}
        struct DefaultObjective {
            std::string namePattern;  // exact match or substring
            bool exactMatch;          // true = exact, false = contains
            int typeIdx;              // 0=SqDev 1=SqOver 2=SqUnder 3=MinDVH 4=MaxDVH
            float doseValue;
            float weight;
            float volumePct;
        };
        
        std::vector<DefaultObjective> defaults = {
            // Lungs: SquaredOverdose at MLD 20 Gy
            {"Poumon_D",          true,  1,  20.0f,  10.0f,  5.0f},
            {"Poumon_G",          true,  1,  20.0f,  10.0f,  5.0f},
            // Heart: MaxDVH V40 < 30%
            {"Coeur",             true,  4,  40.0f,  30.0f, 30.0f},
            // Esophagus: MaxDVH V50 < 30%
            {"Oesophage",         true,  4,  50.0f,  50.0f, 30.0f},
            // Spinal cord: MaxDVH V10Gy <= 0%
            {"Canal_Medullaire",  true,  4,  10.0f, 200.0f,  0.0f},
        };

        for (size_t i = 0; i < ss->getCount(); ++i) {
            const auto* s = ss->getStructure(i);
            if (!s) continue;

            std::string name = s->getName();

            if (s->getType() == "PTV") {
                ObjectiveConfig obj;
                obj.structureName = name;
                obj.typeIdx = 3;        // MinDVH
                obj.doseValue = 66.0f;
                obj.weight = 100.0f;
                obj.volumePct = 98.0f;
                m_objectives.push_back(obj);
            }

            for (const auto& def : defaults) {
                bool match = def.exactMatch 
                    ? (name == def.namePattern)
                    : (name.find(def.namePattern) != std::string::npos);

                if (match) {
                    ObjectiveConfig obj;
                    obj.structureName = name;
                    obj.typeIdx = def.typeIdx;
                    obj.doseValue = def.doseValue;
                    obj.weight = def.weight;
                    obj.volumePct = def.volumePct;
                    m_objectives.push_back(obj);
                }
            }
        }
        m_objectivesInitialized = true;
    }

    // ── Objectives table ── (disabled while optimizing to avoid data races)
    if (m_isOptimizing) ImGui::BeginDisabled();
    if (ImGui::CollapsingHeader("Objectives", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::BeginTable("ObjTable", 6,
                ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
            ImGui::TableSetupColumn("Structure",  ImGuiTableColumnFlags_None, 140.0f);
            ImGui::TableSetupColumn("Type",       ImGuiTableColumnFlags_None, 160.0f);
            ImGui::TableSetupColumn("Dose Gy",  ImGuiTableColumnFlags_None, 80.0f);
            ImGui::TableSetupColumn("Vol %",    ImGuiTableColumnFlags_None, 70.0f);
            ImGui::TableSetupColumn("Weight",     ImGuiTableColumnFlags_None, 80.0f);
            ImGui::TableSetupColumn("##del",      ImGuiTableColumnFlags_None, 25.0f);
            ImGui::TableHeadersRow();

            int toRemove = -1;
            for (size_t i = 0; i < m_objectives.size(); ++i) {
                auto& obj = m_objectives[i];
                ImGui::TableNextRow();
                ImGui::PushID(static_cast<int>(i));

                ImGui::TableSetColumnIndex(0);
                ImGui::SetNextItemWidth(-1);
                if (m_state.patientData && m_state.patientData->getStructureSet()) {
                    const auto* ss = m_state.patientData->getStructureSet();
                    if (ImGui::BeginCombo("##struct", obj.structureName.c_str())) {
                        for (size_t si = 0; si < ss->getCount(); ++si) {
                            const auto* s = ss->getStructure(si);
                            if (!s) continue;
                            bool selected = (obj.structureName == s->getName());
                            if (ImGui::Selectable(s->getName().c_str(), selected)) {
                                obj.structureName = s->getName();
                            }
                            if (selected) ImGui::SetItemDefaultFocus();
                        }
                        ImGui::EndCombo();
                    }
                } else {
                    ImGui::TextUnformatted(obj.structureName.c_str());
                }

                ImGui::TableSetColumnIndex(1);
                ImGui::SetNextItemWidth(-1);
                ImGui::Combo("##type", &obj.typeIdx, kObjectiveTypes, IM_ARRAYSIZE(kObjectiveTypes));

                ImGui::TableSetColumnIndex(2);
                ImGui::SetNextItemWidth(-1);
                ImGui::InputFloat("##dose", &obj.doseValue, 0.0f, 0.0f, "%.1f");

                ImGui::TableSetColumnIndex(3);
                // Volume % only relevant for DVH objectives
                bool isDVH = (obj.typeIdx == 3 || obj.typeIdx == 4);
                if (!isDVH) ImGui::BeginDisabled();
                ImGui::SetNextItemWidth(-1);
                ImGui::InputFloat("##vol", &obj.volumePct, 0.0f, 0.0f, "%.1f");
                if (!isDVH) ImGui::EndDisabled();

                ImGui::TableSetColumnIndex(4);
                ImGui::SetNextItemWidth(-1);
                ImGui::InputFloat("##w", &obj.weight, 0.0f, 0.0f, "%.1f");

                ImGui::TableSetColumnIndex(5);
                if (ImGui::SmallButton("X")) {
                    toRemove = static_cast<int>(i);
                }

                ImGui::PopID();
            }

            if (toRemove >= 0) {
                m_objectives.erase(m_objectives.begin() + toRemove);
            }

            ImGui::EndTable();
        }

        // Add objective button
        if (m_state.patientData && m_state.patientData->getStructureSet()) {
            if (ImGui::Button("+ Add Objective")) {
                ObjectiveConfig obj;
                const auto* ss = m_state.patientData->getStructureSet();
                if (ss->getCount() > 0 && ss->getStructure(0)) {
                    obj.structureName = ss->getStructure(0)->getName();
                }
                m_objectives.push_back(obj);
            }
        }
    }
    if (m_isOptimizing) ImGui::EndDisabled();

    // ── Optimizer settings ──
    ImGui::Spacing();
    if (m_isOptimizing) ImGui::BeginDisabled();
    if (ImGui::CollapsingHeader("Optimizer Settings", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputInt("Max Iterations", &m_maxIterations);
        if (m_maxIterations < 1) m_maxIterations = 1;
        ImGui::InputFloat("Tolerance", &m_tolerance, 1e-7f, 1e-4f, "%.7f");

        ImGui::Separator();
        ImGui::TextUnformatted("Hotspot Control (NTO)");
        ImGui::Checkbox("Enable NTO", &m_ntoEnabled);
        if (!m_ntoEnabled) ImGui::BeginDisabled();
        ImGui::InputFloat("Threshold (%Rx)", &m_ntoThresholdPct, 0.0f, 0.0f, "%.1f");
        ImGui::InputInt("Penalty", &m_ntoPenalty);
        if (!m_ntoEnabled) ImGui::EndDisabled();
    }
    if (m_isOptimizing) ImGui::EndDisabled();

    ImGui::Separator();
    ImGui::Spacing();

    // ── Run Optimization button ──
    if (m_isOptimizing) {
        ImGui::TextColored(getThemeColors().progressText, "Optimizing...");
        if (m_maxIterations > 0) {
            float progress = static_cast<float>(m_currentIteration.load()) / static_cast<float>(m_maxIterations);
            ImGui::ProgressBar(progress, ImVec2(-1, 0),
                ("Iteration " + std::to_string(m_currentIteration.load())).c_str());
        }

        if (m_optimizationDone) {
            m_optThread.join();
            m_isOptimizing = false;
            m_optimizationDone = false;
            m_state.taskRunning = false;

            // Transfer pipeline results to app state
            if (!m_pipelineResult.weights.empty()) {
                m_state.optimizedWeights = std::move(m_pipelineResult.weights);

                if (m_pipelineResult.doseResult) {
                    int optNum = m_state.doseManager.nextOptimizationNumber();
                    m_state.doseManager.addDose(
                        "Optimization #" + std::to_string(optNum),
                        m_pipelineResult.doseResult,
                        m_state.computeGrid);
                    m_state.doseManager.incrementOptimizationCount();
                    m_state.syncSelectedDose();
                    m_state.optimizationJustFinished = true;

                    Logger::info("Max dose: " +
                        std::to_string(m_pipelineResult.doseResult->getMax()) + " Gy");
                }
            }
        }
    } else {
        // Optimization should be possible only when PTV constraints are present to avoid meaningless runs.
        bool canOptimize = m_state.dijComputed() && !m_objectives.empty() &&
                          m_state.patientData && m_state.patientData->getStructureSet() &&
                          std::any_of(m_objectives.begin(), m_objectives.end(),
                                      [this](const ObjectiveConfig& obj) {
                                          const auto* s = m_state.patientData->getStructureSet()
                                              ->getStructureByName(obj.structureName);
                                          return s && s->isTarget();
                                      });
        if (!canOptimize) {
            ImGui::TextColored(getThemeColors().failText, "No valid PTV constraints.");
            ImGui::BeginDisabled();
        }

        if (ImGui::Button("Run Optimization", ImVec2(-1, 30))) {
            // Join any previous thread before starting a new one
            if (m_optThread.joinable()) m_optThread.join();

            m_state.resetOptimization();
            m_isOptimizing = true;
            m_optimizationDone = false;
            m_state.taskRunning = true;
            m_currentIteration = 0;

            // Snapshot objectives to avoid data race with render thread
            auto objectivesCopy = m_objectives;
            int maxIter = m_maxIterations;
            float tolerance = m_tolerance;
            bool ntoEnabled = m_ntoEnabled;
            float ntoThresholdPct = m_ntoThresholdPct;
            int ntoPenalty = m_ntoPenalty;

            m_optThread = std::thread([this, objectivesCopy, maxIter, tolerance,
                                       ntoEnabled, ntoThresholdPct, ntoPenalty]() {
                try {
                    const auto& ctGrid = m_state.patientData->getGrid();
                    const auto& doseGrid = *m_state.computeGrid;

                    // Convert GUI objectives to ObjectiveSpecs
                    std::vector<ObjectiveSpec> specs;
                    for (const auto& cfg : objectivesCopy) {
                        ObjectiveSpec spec;
                        spec.structurePattern = cfg.structureName;
                        spec.exactMatch = true;
                        spec.typeIdx = cfg.typeIdx;
                        spec.doseValue = static_cast<double>(cfg.doseValue);
                        spec.weight = static_cast<double>(cfg.weight);
                        spec.volumePct = static_cast<double>(cfg.volumePct);
                        specs.push_back(spec);
                    }

                    auto objectives = ObjectiveBuilder::buildFromSpecs(
                        specs, *m_state.patientData, ctGrid, doseGrid);

                    // Detect prescription dose from target objectives
                    double prescriptionDose = 0.0;
                    const auto* ss = m_state.patientData->getStructureSet();
                    for (const auto& cfg : objectivesCopy) {
                        const auto* structure = ss->getStructureByName(cfg.structureName);
                        if (structure && structure->isTarget() &&
                            cfg.doseValue > prescriptionDose) {
                            prescriptionDose = cfg.doseValue;
                        }
                    }

                    OptimizationConfig config;
                    config.maxIterations = maxIter;
                    config.tolerance = static_cast<double>(tolerance);
                    config.targetDose = prescriptionDose > 0 ? prescriptionDose : 66.0;
                    config.ntoEnabled = ntoEnabled;
                    config.ntoThresholdPct = static_cast<double>(ntoThresholdPct) / 100.0;
                    config.ntoPenalty = static_cast<double>(ntoPenalty);

                    m_pipelineResult = OptimizationPipeline::runWithObjectives(
                        *m_state.dij, config, std::move(objectives),
                        *m_state.patientData, doseGrid);

                    m_optStatusMessage = (m_pipelineResult.converged ? "Converged" : "Max iter reached") +
                        std::string(" in ") + std::to_string(m_pipelineResult.iterations) +
                        " iterations | Obj=" + std::to_string(m_pipelineResult.finalObjective);
                } catch (const std::exception& e) {
                    Logger::error("Optimization failed: " + std::string(e.what()));
                    m_optStatusMessage = "Error: " + std::string(e.what());
                }

                m_optimizationDone = true;
            });
        }

        if (!canOptimize) ImGui::EndDisabled();
    }

    // ── Status ──
    if (m_state.optimizationDone()) {
        ImGui::Spacing();
        ImGui::TextColored(getThemeColors().passText, "%s", m_optStatusMessage.c_str());
        if (m_state.doseResult) {
            ImGui::Text("Max dose: %.2f Gy | Mean dose: %.2f Gy",
                m_state.doseResult->getMax(), m_state.doseResult->getMean());
        }
    }

    ImGui::End();
}

} // namespace optirad
