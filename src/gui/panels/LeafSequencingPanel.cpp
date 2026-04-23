#include "LeafSequencingPanel.hpp"
#include "Theme.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>

namespace optirad {

LeafSequencingPanel::LeafSequencingPanel(GuiAppState& state) : m_state(state) {}

LeafSequencingPanel::~LeafSequencingPanel() {
    if (m_thread.joinable()) m_thread.join();
}

void LeafSequencingPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Leaf Sequencing", &m_visible);

    if (!m_state.optimizationDone()) {
        ImGui::TextDisabled("Run optimization first to enable leaf sequencing.");
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
                        *m_state.dij,
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
