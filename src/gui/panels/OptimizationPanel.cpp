#include "OptimizationPanel.hpp"
#include "optimization/OptimizerFactory.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/objectives/SquaredDeviation.hpp"
#include "optimization/objectives/SquaredOverdose.hpp"
#include "optimization/objectives/SquaredUnderdose.hpp"
#include "optimization/objectives/DVHObjective.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "dose/PlanAnalysis.hpp"
#include "geometry/Grid.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <chrono>
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

    // ── Initialize objectives from structures ──
    if (!m_objectivesInitialized && m_state.patientData && m_state.patientData->getStructureSet()) {
        const auto* ss = m_state.patientData->getStructureSet();
        for (size_t i = 0; i < ss->getCount(); ++i) {
            const auto* s = ss->getStructure(i);
            if (!s) continue;
            ObjectiveConfig obj;
            obj.structureName = s->getName();
            std::string type = s->getType();
            std::string name = s->getName();
            bool isTarget = (type.find("TARGET") != std::string::npos ||
                           type.find("PTV") != std::string::npos ||
                           type.find("CTV") != std::string::npos ||
                           type.find("GTV") != std::string::npos ||
                           name.find("PTV") != std::string::npos ||
                           name.find("CTV") != std::string::npos ||
                           name.find("GTV") != std::string::npos);
            if (isTarget) {
                obj.typeIdx = 0; // SquaredDeviation
                obj.doseValue = 60.0f;
                obj.weight = 100.0f;
                obj.volumePct = 95.0f;
            } else {
                obj.typeIdx = 1; // SquaredOverdose (OAR)
                obj.doseValue = 30.0f;
                obj.weight = 1.0f;
                obj.volumePct = 5.0f;
            }
            m_objectives.push_back(obj);
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
            ImGui::TableSetupColumn("Dose (Gy)",  ImGuiTableColumnFlags_None, 80.0f);
            ImGui::TableSetupColumn("Vol (%)",    ImGuiTableColumnFlags_None, 70.0f);
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
                ImGui::DragFloat("##dose", &obj.doseValue, 0.5f, 0.0f, 120.0f, "%.1f");

                ImGui::TableSetColumnIndex(3);
                // Volume % only relevant for DVH objectives
                bool isDVH = (obj.typeIdx == 3 || obj.typeIdx == 4);
                if (!isDVH) ImGui::BeginDisabled();
                ImGui::SetNextItemWidth(-1);
                ImGui::DragFloat("##vol", &obj.volumePct, 1.0f, 0.0f, 100.0f, "%.0f");
                if (!isDVH) ImGui::EndDisabled();

                ImGui::TableSetColumnIndex(4);
                ImGui::SetNextItemWidth(-1);
                ImGui::DragFloat("##w", &obj.weight, 0.1f, 0.0f, 1000.0f, "%.1f");

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
    if (ImGui::CollapsingHeader("Optimizer Settings")) {
        ImGui::InputInt("Max Iterations", &m_maxIterations);
        if (m_maxIterations < 1) m_maxIterations = 1;
        ImGui::InputFloat("Tolerance", &m_tolerance, 1e-7f, 1e-4f, "%.7f");

        ImGui::Separator();
        ImGui::TextUnformatted("Hotspot Control (NTO)");
        ImGui::Checkbox("Enable NTO", &m_ntoEnabled);
        if (!m_ntoEnabled) ImGui::BeginDisabled();
        ImGui::DragFloat("Threshold (%Rx)", &m_ntoThresholdPct, 0.5f, 90.0f, 120.0f, "%.0f%%");
        ImGui::DragFloat("Penalty", &m_ntoPenalty, 100.0f, 0.0f, 100000.0f, "%.0f");
        if (!m_ntoEnabled) ImGui::EndDisabled();
    }
    if (m_isOptimizing) ImGui::EndDisabled();

    ImGui::Separator();
    ImGui::Spacing();

    // ── Run Optimization button ──
    if (m_isOptimizing) {
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Optimizing...");
        if (m_maxIterations > 0) {
            float progress = static_cast<float>(m_currentIteration.load()) / static_cast<float>(m_maxIterations);
            ImGui::ProgressBar(progress, ImVec2(-1, 0),
                ("Iteration " + std::to_string(m_currentIteration.load())).c_str());
        }

        if (m_optimizationDone) {
            m_optThread.join();
            m_isOptimizing = false;
            m_optimizationDone = false;

            // Compute forward dose
            if (!m_state.optimizedWeights.empty()) {
                try {
                    auto engine = DoseEngineFactory::create("PencilBeam");
                    auto dose = engine->calculateDose(*m_state.dij, m_state.optimizedWeights, *m_state.doseGrid);
                    m_state.doseResult = std::make_shared<DoseMatrix>(std::move(dose));
                    Logger::info("Forward dose computed. Max dose: " +
                        std::to_string(m_state.doseResult->getMax()) + " Gy");

                    // Post-optimization plan analysis
                    auto stats = PlanAnalysis::computeStats(
                        *m_state.doseResult, *m_state.patientData, *m_state.doseGrid);
                    std::ostringstream oss;
                    PlanAnalysis::print(stats, oss);
                    Logger::info(oss.str());
                } catch (const std::exception& e) {
                    Logger::error("Forward dose failed: " + std::string(e.what()));
                }
            }
        }
    } else {
        bool canOptimize = m_state.dijComputed() && !m_objectives.empty();
        if (!canOptimize) ImGui::BeginDisabled();

        if (ImGui::Button("Run Optimization", ImVec2(-1, 30))) {
            // Join any previous thread before starting a new one
            if (m_optThread.joinable()) m_optThread.join();

            m_state.resetOptimization();
            m_isOptimizing = true;
            m_optimizationDone = false;
            m_currentIteration = 0;

            // Snapshot objectives to avoid data race with render thread
            auto objectivesCopy = m_objectives;
            int maxIter = m_maxIterations;
            float tolerance = m_tolerance;
            bool ntoEnabled = m_ntoEnabled;
            float ntoThresholdPct = m_ntoThresholdPct;
            float ntoPenalty = m_ntoPenalty;

            m_optThread = std::thread([this, objectivesCopy, maxIter, tolerance,
                                       ntoEnabled, ntoThresholdPct, ntoPenalty]() {
                auto start = std::chrono::steady_clock::now();

                try {
                    const auto* ss = m_state.patientData->getStructureSet();
                    const auto& ctGrid = m_state.patientData->getGrid();
                    const auto& doseGrid = *m_state.doseGrid;

                    // ── Pre-optimization logging ──
                    Logger::info("=== Pre-Optimization Summary ===");
                    Logger::info("Structures: " + std::to_string(ss->getCount()) +
                                " | Bixels: " + std::to_string(m_state.dij->getNumBixels()) +
                                " | Dose voxels: " + std::to_string(m_state.dij->getNumVoxels()) +
                                " | NNZ: " + std::to_string(m_state.dij->getNumNonZeros()));

                    // Build objectives
                    std::vector<std::unique_ptr<ObjectiveFunction>> ownedObjs;
                    std::vector<ObjectiveFunction*> objPtrs;

                    for (const auto& cfg : objectivesCopy) {
                        const auto* structure = ss->getStructureByName(cfg.structureName);
                        if (!structure) continue;

                        auto doseIndices = Grid::mapVoxelIndices(
                            ctGrid, doseGrid, structure->getVoxelIndices());

                        std::unique_ptr<ObjectiveFunction> obj;
                        std::string objDesc;

                        switch (cfg.typeIdx) {
                            case 0: {
                                auto sd = std::make_unique<SquaredDeviation>();
                                sd->setPrescribedDose(static_cast<double>(cfg.doseValue));
                                objDesc = "SquaredDeviation @ " + std::to_string(cfg.doseValue) + " Gy";
                                obj = std::move(sd);
                                break;
                            }
                            case 1: {
                                auto so = std::make_unique<SquaredOverdose>();
                                so->setMaxDose(static_cast<double>(cfg.doseValue));
                                objDesc = "SquaredOverdose @ " + std::to_string(cfg.doseValue) + " Gy";
                                obj = std::move(so);
                                break;
                            }
                            case 2: {
                                auto su = std::make_unique<SquaredUnderdose>();
                                su->setMinDose(static_cast<double>(cfg.doseValue));
                                objDesc = "SquaredUnderdose @ " + std::to_string(cfg.doseValue) + " Gy";
                                obj = std::move(su);
                                break;
                            }
                            case 3: {
                                auto dvh = std::make_unique<DVHObjective>();
                                dvh->setType(DVHObjective::Type::MIN_DVH);
                                dvh->setDoseThreshold(static_cast<double>(cfg.doseValue));
                                dvh->setVolumeFraction(static_cast<double>(cfg.volumePct) / 100.0);
                                objDesc = "MinDVH D" + std::to_string(static_cast<int>(cfg.volumePct)) +
                                         "% >= " + std::to_string(cfg.doseValue) + " Gy";
                                obj = std::move(dvh);
                                break;
                            }
                            case 4: {
                                auto dvh = std::make_unique<DVHObjective>();
                                dvh->setType(DVHObjective::Type::MAX_DVH);
                                dvh->setDoseThreshold(static_cast<double>(cfg.doseValue));
                                dvh->setVolumeFraction(static_cast<double>(cfg.volumePct) / 100.0);
                                objDesc = "MaxDVH V" + std::to_string(cfg.doseValue) +
                                         "Gy <= " + std::to_string(static_cast<int>(cfg.volumePct)) + "%";
                                obj = std::move(dvh);
                                break;
                            }
                        }

                        if (obj) {
                            obj->setWeight(static_cast<double>(cfg.weight));
                            obj->setStructure(structure);
                            obj->setVoxelIndices(doseIndices);

                            Logger::info("  " + cfg.structureName + " (" + structure->getType() +
                                        ") -> " + objDesc + " | w=" +
                                        std::to_string(cfg.weight) +
                                        " | " + std::to_string(doseIndices.size()) + " dose voxels");

                            objPtrs.push_back(obj.get());
                            ownedObjs.push_back(std::move(obj));
                        }
                    }

                    Logger::info("Total objectives: " + std::to_string(objPtrs.size()));
                    Logger::info("Max iterations: " + std::to_string(maxIter) +
                                " | Tolerance: " + std::to_string(tolerance));

                    // Detect prescription dose for NTO from ALL objective types
                    // Use the max dose value from any target structure objective
                    double prescriptionDose = 0.0;
                    for (const auto& cfg : objectivesCopy) {
                        // Check if this is a target structure
                        const auto* structure = ss->getStructureByName(cfg.structureName);
                        if (!structure) continue;
                        std::string sType = structure->getType();
                        std::string sName = structure->getName();
                        bool isTarget = (sType.find("TARGET") != std::string::npos ||
                                        sType.find("PTV") != std::string::npos ||
                                        sType.find("CTV") != std::string::npos ||
                                        sType.find("GTV") != std::string::npos ||
                                        sName.find("PTV") != std::string::npos ||
                                        sName.find("CTV") != std::string::npos ||
                                        sName.find("GTV") != std::string::npos);
                        if (isTarget && cfg.doseValue > prescriptionDose) {
                            prescriptionDose = cfg.doseValue;
                        }
                    }

                    // Create optimizer
                    auto optimizer = OptimizerFactory::create("LBFGS");
                    optimizer->setMaxIterations(maxIter);
                    optimizer->setTolerance(static_cast<double>(tolerance));

                    // Configure NTO (hotspot control) if enabled and we have a prescription dose
                    auto* lbfgs = dynamic_cast<LBFGSOptimizer*>(optimizer.get());
                    if (lbfgs && ntoEnabled && prescriptionDose > 0) {
                        lbfgs->setPrescriptionDose(prescriptionDose);
                        lbfgs->setHotspotThreshold(static_cast<double>(ntoThresholdPct) / 100.0);
                        lbfgs->setHotspotPenalty(static_cast<double>(ntoPenalty));
                        Logger::info("NTO enabled: Rx=" + std::to_string(prescriptionDose) +
                                    " Gy, threshold=" + std::to_string(ntoThresholdPct) +
                                    "%, penalty=" + std::to_string(ntoPenalty));
                    }

                    // Run
                    std::vector<Constraint> constraints;
                    auto result = optimizer->optimize(*m_state.dij, objPtrs, constraints);

                    m_state.optimizedWeights = std::move(result.weights);

                    auto end = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration<double>(end - start).count();
                    m_optStatusMessage = (result.converged ? "Converged" : "Max iter reached") +
                        std::string(" in ") + std::to_string(result.iterations) +
                        " iterations (" + std::to_string(elapsed).substr(0, 5) + "s)" +
                        " | Obj=" + std::to_string(result.finalObjective);

                    Logger::info("Optimization: " + m_optStatusMessage);
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
        ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "%s", m_optStatusMessage.c_str());
        if (m_state.doseResult) {
            ImGui::Text("Max dose: %.2f Gy | Mean dose: %.2f Gy",
                m_state.doseResult->getMax(), m_state.doseResult->getMean());
        }
    }

    ImGui::End();
}

} // namespace optirad
