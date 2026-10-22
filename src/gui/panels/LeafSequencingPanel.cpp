#include "LeafSequencingPanel.hpp"
#include "Theme.hpp"
#include "optimization/ObjectiveBuilder.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <algorithm>
#include <cmath>

namespace optirad {

LeafSequencingPanel::LeafSequencingPanel(GuiAppState& state) : m_state(state) {}

LeafSequencingPanel::~LeafSequencingPanel() {
    if (m_thread.joinable()) m_thread.join();
    if (m_arcThread.joinable()) m_arcThread.join();
}

void LeafSequencingPanel::renderArc() {
    const bool running = m_arcRunning.load();

    ImGui::TextUnformatted("VMAT: arc sequencing and aperture optimization");
    ImGui::TextDisabled("DAO: apertures are optimized directly against the objectives at increasing resolution.\n"
                        "Sweep: the optimized fluence is delivered with sliding windows, then refined.");

    if (running) ImGui::BeginDisabled();
    if (ImGui::CollapsingHeader("Arc Sequencer Settings", ImGuiTreeNodeFlags_DefaultOpen)) {
        static const char* kStrategies[] = {"Fluence sweep + refinement", "Progressive-resolution DAO (Eclipse-like)"};
        ImGui::Combo("Strategy", &m_arcStrategy, kStrategies, IM_ARRAYSIZE(kStrategies));
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("DAO: conformal apertures at a few control points are optimized directly against the\n"
                              "objectives (leaf positions and MU), then control points are inserted and re-optimized\n"
                              "until the plan's control point spacing. Few MU, smooth leaf motion, no gantry slowdown.");
        const bool dao = m_arcStrategy == static_cast<int>(ArcStrategy::ProgressiveDao);

        if (dao) {
            ImGui::InputInt("Iterations per resolution", &m_arcDaoIterations);
            m_arcDaoIterations = std::clamp(m_arcDaoIterations, 10, 5000);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("The last (finest) resolution gets twice as many.");
            ImGui::InputFloat("First spacing (deg)", &m_arcDaoCoarsest, 0.0f, 0.0f, "%.1f");
            m_arcDaoCoarsest = std::max(m_arcDaoCoarsest, 1.0f);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Control point spacing of the first resolution; halved until the plan's spacing.");
            ImGui::InputFloat("MU penalty", &m_arcDaoMuPenalty, 0.05f, 0.25f, "%.2f");
            m_arcDaoMuPenalty = std::max(m_arcDaoMuPenalty, 0.0f);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Cost of the total MU as a fraction of the objective at the second resolution. 0 = none.");
        } else {
            static const char* kMethods[] = {"Auto (sweep when possible)", "Sweep (sliding window)",
                                             "Intensity levels", "Conformal (seed only)"};
            ImGui::Combo("Sequencing", &m_arcMethod, kMethods, IM_ARRAYSIZE(kMethods));
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Sweep needs an interdigitating MLC and at least 2 control points per sector.");
            if (m_arcMethod == static_cast<int>(ArcSequencerOptions::Method::LevelSet)) {
                ImGui::SliderInt("Max intensity levels", &m_arcLevels, 3, 30);
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Upper bound; each sector is limited to its number of control points.");
            } else {
                ImGui::InputFloat("Max fluence rise", &m_arcRiseRatio, 0.05f, 0.25f, "%.2f");
                if (m_arcRiseRatio < 0.0f) m_arcRiseRatio = 0.0f;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Rows whose upward fluence steps exceed this multiple of the sector's peak are smoothed.\n"
                                      "Lower = fewer MU, less modulation. 0 = off.");
                ImGui::InputFloat("Max sweep fraction", &m_arcSweepFraction, 0.05f, 0.1f, "%.2f");
                m_arcSweepFraction = std::clamp(m_arcSweepFraction, 0.05f, 0.95f);
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Largest share of a sector's time the leaves may spend sweeping across the field.");
            }
        }

        ImGui::Checkbox("Integrate leaf motion", &m_arcIntegrateMotion);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Dose of a control point = average over the leaf motion in its slice.");
        ImGui::InputFloat("Hot-spot threshold (%Rx)", dao ? &m_arcDaoHotspotThresholdPct : &m_arcHotspotThresholdPct,
                          0.0f, 0.0f, "%.1f");
        ImGui::InputInt("Hot-spot penalty", dao ? &m_arcDaoHotspotPenalty : &m_arcHotspotPenalty);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Stricter than the fluence stage: delivery errors add hot spots to the arc dose.\n"
                              "DAO tightens from the fluence stage's setting to this one over its resolutions.");
        ImGui::InputFloat("MU smoothing", &m_arcMuSmoothing, 0.01f, 0.1f, "%.3f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Penalty on MU differences between neighbouring control points.");

        if (!dao) {
            ImGui::Checkbox("Optimize aperture weights", &m_arcOptimizeWeights);
            ImGui::Checkbox("Optimize leaf positions", &m_arcOptimizeLeaves);
            if (m_arcOptimizeLeaves) {
                ImGui::InputInt("Leaf iterations", &m_arcLeafIterations);
                m_arcLeafIterations = std::clamp(m_arcLeafIterations, 1, 2000);
            }
            if (!m_arcOptimizeWeights) ImGui::BeginDisabled();
            ImGui::InputInt("Weight iterations", &m_arcIterations);
            if (m_arcIterations < 1) m_arcIterations = 1;
            ImGui::InputFloat("MU penalty", &m_arcMuPenalty, 0.05f, 0.25f, "%.2f");
            if (m_arcMuPenalty < 0.0f) m_arcMuPenalty = 0.0f;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Cost of the total MU as a fraction of the initial objective. 0 = none.");
            ImGui::InputFloat("Max MU growth", &m_arcMaxWeightFactor, 0.5f, 1.0f, "%.1f");
            if (m_arcMaxWeightFactor < 1.0f) m_arcMaxWeightFactor = 1.0f;
            if (!m_arcOptimizeWeights) ImGui::EndDisabled();
        }
    }
    if (running) ImGui::EndDisabled();

    ImGui::Separator();
    ImGui::Spacing();

    if (running) {
        {
            std::lock_guard<std::mutex> lock(m_arcMutex);
            ImGui::TextColored(getThemeColors().progressText, "%s...", m_arcStage.c_str());
        }
        ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));

        if (m_arcDone) {
            m_arcThread.join();
            m_arcRunning = false;
            m_arcDone = false;
            m_state.taskRunning = false;

            if (m_arcPending) {
                const int seqNum = m_state.nextLeafSeqNumber();
                const int doseId = m_state.doseStore.addDose(
                    "VMAT Arc #" + std::to_string(seqNum), m_arcPending->dose, m_state.computeGrid);
                m_state.incrementLeafSeqCount();
                m_state.arcCache[doseId] = {m_arcPending, m_state.activeOptDoseId};
                m_state.syncSelectedDose();
                m_arcPending.reset();
            }
        }
    } else {
        const bool hasSpecs = !m_state.optimizationSpecs.empty() && m_state.computeGrid;
        if (!hasSpecs) ImGui::BeginDisabled();
        if (ImGui::Button("Run Arc Sequencing", ImVec2(-1, 30))) {
            if (m_arcThread.joinable()) m_arcThread.join();
            m_state.resetLeafSequence();
            m_arcRunning = true;
            m_arcDone = false;
            m_state.taskRunning = true;
            m_arcPending.reset();
            m_arcStatus.clear();

            ArcOptimizationConfig settings;
            settings.strategy = static_cast<ArcStrategy>(m_arcStrategy);
            settings.daoIterations = m_arcDaoIterations;
            settings.daoCoarsestSpacingDeg = static_cast<double>(m_arcDaoCoarsest);
            settings.daoHotspotThresholdPct = 0.01 * static_cast<double>(m_arcDaoHotspotThresholdPct);
            settings.daoHotspotPenalty = static_cast<double>(m_arcDaoHotspotPenalty);
            settings.daoMuPenalty = static_cast<double>(m_arcDaoMuPenalty);
            settings.sequencer.method = static_cast<ArcSequencerOptions::Method>(m_arcMethod);
            settings.sequencer.numLevels = m_arcLevels;
            settings.sequencer.maxRiseRatio = static_cast<double>(m_arcRiseRatio);
            settings.sequencer.maxSweepFraction = static_cast<double>(m_arcSweepFraction);
            settings.integrateLeafMotion = m_arcIntegrateMotion;
            settings.optimizeWeights = m_arcOptimizeWeights;
            settings.optimizeLeaves = m_arcOptimizeLeaves;
            settings.leafIterations = m_arcLeafIterations;
            settings.weightIterations = m_arcIterations;
            settings.muSmoothingWeight = static_cast<double>(m_arcMuSmoothing);
            settings.hotspotThresholdPct = 0.01 * static_cast<double>(m_arcHotspotThresholdPct);
            settings.hotspotPenalty = static_cast<double>(m_arcHotspotPenalty);
            settings.muPenalty = static_cast<double>(m_arcMuPenalty);
            settings.maxWeightFactor = static_cast<double>(m_arcMaxWeightFactor);

            m_arcThread = std::thread([this, settings]() {
                try {
                    const Grid& doseGrid = *m_state.computeGrid;
                    ArcOptimizationConfig cfg = settings;
                    cfg.base = m_state.optimizationConfig;
                    const double res = m_state.plan->getMachine().getGeometry().leafPositionResolution;
                    if (res > 0.0) cfg.sequencer.leafPositionResolution = res;
                    cfg.computeBackend = m_state.computeBackend;
                    cfg.gpuMemoryBudgetMiB = m_state.gpuMemoryBudgetMiB;
                    if (m_state.doseOp) cfg.engineOptions = m_state.doseOp->engineOptions();

                    auto objectives = ObjectiveBuilder::buildFromSpecs(
                        m_state.optimizationSpecs, *m_state.patientData, m_state.patientData->getGrid(), doseGrid);

                    m_arcPending = std::make_shared<ArcPipelineResult>(ArcPipeline::run(
                        m_state.optimizedWeights, *m_state.stf, *m_state.plan, *m_state.patientData, doseGrid,
                        std::move(objectives), cfg,
                        [this](const IterationInfo& info) { m_arcIter.store(info.iteration); },
                        [this](const std::string& stage) {
                            std::lock_guard<std::mutex> lock(m_arcMutex);
                            m_arcStage = stage;
                        }));

                    const auto& seq = m_arcPending->sequence;
                    char buf[256];
                    snprintf(buf, sizeof(buf), "%zu control points | %.0f MU | %.0f s | %d slowed | %.1f mm leaf travel/cp",
                             seq.controlPoints.size(), seq.totalMU, seq.totalTimeSec, seq.slowedControlPoints,
                             seq.meanLeafTravel);
                    m_arcStatus = buf;
                } catch (const std::exception& e) {
                    Logger::error("Arc sequencing failed: " + std::string(e.what()));
                    m_arcStatus = "Error: " + std::string(e.what());
                    m_arcPending.reset();
                }
                m_arcDone = true;
            });
        }
        if (!hasSpecs) {
            ImGui::EndDisabled();
            ImGui::TextDisabled("Run the fluence optimization in this session first (its objectives are reused).");
        }
    }

    // ── Results ──
    if (!m_arcStatus.empty() && !running) {
        ImGui::Spacing();
        ImGui::TextColored(m_state.arcResult ? getThemeColors().passText : getThemeColors().failText,
                           "%s", m_arcStatus.c_str());
    }
    if (m_state.arcResult && !running) {
        const auto& seq = m_state.arcResult->sequence;
        ImGui::Text("Leaf refinement: %d steps | %.3e -> %.3e",
            m_state.arcResult->leafIterations, m_state.arcResult->initialLeafObjective,
            m_state.arcResult->finalLeafObjective);
        ImGui::Text("Slowest gantry speed: %.2f deg/s | max leaf speed: %.1f mm/s",
                    seq.minGantrySpeed, seq.maxLeafSpeed);
        if (m_arcOptimizeWeights)
            ImGui::Text("Aperture optimization: %d iterations, objective %.3e",
                        m_state.arcResult->iterations, m_state.arcResult->finalObjective);

        if (ImGui::CollapsingHeader("Control Points")) {
            if (ImGui::BeginTable("CpTable", 6,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit |
                    ImGuiTableFlags_ScrollY, ImVec2(0, 220))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("#",        ImGuiTableColumnFlags_None, 35.0f);
                ImGui::TableSetupColumn("Gantry",   ImGuiTableColumnFlags_None, 60.0f);
                ImGui::TableSetupColumn("MU",       ImGuiTableColumnFlags_None, 60.0f);
                ImGui::TableSetupColumn("Time s",   ImGuiTableColumnFlags_None, 60.0f);
                ImGui::TableSetupColumn("deg/s",    ImGuiTableColumnFlags_None, 55.0f);
                ImGui::TableSetupColumn("Leaf mm/s", ImGuiTableColumnFlags_None, 70.0f);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < seq.controlPoints.size(); ++i) {
                    const auto& cp = seq.controlPoints[i];
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%zu", i);
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%.1f", cp.geometry.gantryAngle);
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.2f", cp.mu);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%.2f", cp.timeSec);
                    ImGui::TableSetColumnIndex(4); ImGui::Text("%.2f", cp.gantrySpeed);
                    ImGui::TableSetColumnIndex(5); ImGui::Text("%.1f", cp.leafSpeed);
                }
                ImGui::EndTable();
            }
        }
    }
}

void LeafSequencingPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Leaf Sequencing", &m_visible);

    if (!m_state.optimizationDone()) {
        ImGui::TextDisabled("Run optimization first to enable leaf sequencing.");
        ImGui::End();
        return;
    }

    if (m_state.plan && m_state.plan->isVmat()) {
        renderArc();
        ImGui::End();
        return;
    }

    // Auto-restore display state when selected dose changes
    if (!m_isRunning) {
        int ver = m_state.doseStore.version();
        if (ver != m_doseVersion) {
            m_doseVersion = ver;
            // Find the active seq cache entry for this dose selection
            // Check if selected dose is a deliverable dose directly
            auto* sel = m_state.doseStore.getSelected();
            int selId = sel ? sel->id : -1;
            auto it = m_displayCache.find(selId);
            if (it != m_displayCache.end()) {
                m_statusMessage = it->second.statusMessage;
            } else {
                // Check if an optimization dose is selected that has a linked seq
                bool found = false;
                for (const auto& [seqDoseId, seqEntry] : m_state.seqCache) {
                    if (seqEntry.linkedOptDoseId == selId) {
                        auto dit = m_displayCache.find(seqDoseId);
                        if (dit != m_displayCache.end()) {
                            m_statusMessage = dit->second.statusMessage;
                            found = true;
                        }
                        break;
                    }
                }
                if (!found) {
                    m_statusMessage.clear();
                }
            }
        }
    }

    // ── Settings ──
    if (m_isRunning) ImGui::BeginDisabled();
    if (ImGui::CollapsingHeader("Sequencer Settings", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SliderInt("Intensity Levels", &m_numLevels, 5, 30);
        ImGui::InputFloat("Min Segment MU", &m_minSegmentMU, 0.1f, 1.0f, "%.1f");
        if (m_minSegmentMU < 0.0f) m_minSegmentMU = 0.0f;
    }
    if (m_isRunning) ImGui::EndDisabled();

    ImGui::Separator();
    ImGui::Spacing();

    // ── Run button / progress ──
    if (m_isRunning) {
        ImGui::TextColored(getThemeColors().progressText, "Sequencing...");
        int cur = m_currentBeam.load();
        int tot = m_totalBeams.load();
        if (tot > 0) {
            float progress = static_cast<float>(cur) / static_cast<float>(tot);
            ImGui::ProgressBar(progress, ImVec2(-1, 0),
                ("Beam " + std::to_string(cur) + " / " + std::to_string(tot)).c_str());
        }

        if (m_isDone) {
            m_thread.join();
            m_isRunning = false;
            m_isDone = false;
            m_state.taskRunning = false;

            if (!m_pipelineResult.beamSequences.empty()) {
                m_state.leafSequences = std::move(m_pipelineResult.beamSequences);
                m_state.deliverableWeights = std::move(m_pipelineResult.deliverableWeights);
                m_state.deliverableStats = std::move(m_pipelineResult.deliverableStats);

                // Add deliverable dose to DoseManager
                if (m_pipelineResult.deliverableDose) {
                    m_state.deliverableDoseResult = m_pipelineResult.deliverableDose;
                    int seqNum = m_state.nextLeafSeqNumber();
                    int doseId = m_state.doseStore.addDose(
                        "Leaf Sequencing #" + std::to_string(seqNum),
                        m_pipelineResult.deliverableDose,
                        m_state.computeGrid);
                    m_state.incrementLeafSeqCount();

                    // Cache leaf seq results BEFORE syncSelectedDose
                    m_state.cacheLeafSequencing(doseId, m_state.activeOptDoseId);

                    // Cache display state
                    m_displayCache[doseId] = {m_statusMessage};

                    m_state.syncSelectedDose();
                }
            }
        }
    } else {
        if (ImGui::Button("Run Leaf Sequencing", ImVec2(-1, 30))) {
            if (m_thread.joinable()) m_thread.join();

            m_state.resetLeafSequence();
            m_isRunning = true;
            m_isDone = false;
            m_state.taskRunning = true;
            m_currentBeam = 0;
            m_totalBeams = 0;

            int numLevels = m_numLevels;
            float minMU = m_minSegmentMU;

            m_thread = std::thread([this, numLevels, minMU]() {
                try {
                    LeafSequencerOptions opts;
                    opts.numLevels = numLevels;
                    opts.minSegmentMU = static_cast<double>(minMU);

                    // Use leaf position resolution from machine geometry if available
                    if (m_state.plan) {
                        const auto& mlc = m_state.plan->getMachine().getGeometry();
                        if (mlc.leafPositionResolution > 0.0) {
                            opts.leafPositionResolution = mlc.leafPositionResolution;
                        }
                    }

                    m_pipelineResult = LeafSequencingPipeline::run(
                        m_state.optimizedWeights,
                        *m_state.stf,
                        *m_state.doseOp,
                        *m_state.plan,
                        *m_state.patientData,
                        *m_state.displayGrid,
                        opts,
                        [this](int cur, int tot) {
                            m_currentBeam.store(cur);
                            m_totalBeams.store(tot);
                        });

                    m_statusMessage = std::to_string(m_pipelineResult.totalSegments) +
                        " segments | " + std::to_string(static_cast<int>(m_pipelineResult.totalMU)) +
                        " MU | fidelity " +
                        std::to_string(m_pipelineResult.meanFidelity).substr(0, 6);
                } catch (const std::exception& e) {
                    Logger::error("Leaf sequencing failed: " + std::string(e.what()));
                    m_statusMessage = "Error: " + std::string(e.what());
                }
                m_isDone = true;
            });
        }
    }

    // ── Results ──
    if (m_state.leafSequenceDone()) {
        ImGui::Spacing();
        ImGui::TextColored(getThemeColors().passText, "%s", m_statusMessage.c_str());

        if (ImGui::CollapsingHeader("Per-Beam Results", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::BeginTable("SeqTable", 4,
                    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupColumn("Beam",     ImGuiTableColumnFlags_None, 50.0f);
                ImGui::TableSetupColumn("Segments", ImGuiTableColumnFlags_None, 70.0f);
                ImGui::TableSetupColumn("MU",       ImGuiTableColumnFlags_None, 80.0f);
                ImGui::TableSetupColumn("Fidelity", ImGuiTableColumnFlags_None, 80.0f);
                ImGui::TableHeadersRow();

                for (size_t i = 0; i < m_state.leafSequences.size(); ++i) {
                    const auto& seq = m_state.leafSequences[i];
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    ImGui::Text("%zu", i);
                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%zu", seq.segments.size());
                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.1f", seq.totalMU);
                    ImGui::TableSetColumnIndex(3);
                    ImGui::Text("%.4f", seq.fluenceFidelity);
                }

                ImGui::EndTable();
            }
        }
    }

    ImGui::End();
}

} // namespace optirad
