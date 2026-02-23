#include "PhaseSpacePanel.hpp"
#include "../views/renderers/PhaseSpaceRenderer.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "phsp/PhaseSpaceData.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <chrono>
#include <cmath>
#include <numeric>

namespace optirad {

PhaseSpacePanel::PhaseSpacePanel(GuiAppState& state) : m_state(state) {}

void PhaseSpacePanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Phase Space", &m_visible);

    // Gate: need a phase-space plan
    if (!m_state.planCreated()) {
        ImGui::TextDisabled("Create a plan with a phase-space machine first.");
        ImGui::End();
        return;
    }

    if (!m_state.isPhaseSpaceMachine()) {
        ImGui::TextDisabled("Current machine is Generic (pencil-beam).");
        ImGui::TextDisabled("Phase-space panel is for machines with IAEA PSF data.");
        ImGui::End();
        return;
    }

    // Show machine info
    const auto& machine = m_state.plan->getMachine();
    const auto& geom = machine.getGeometry();
    ImGui::Text("Machine: %s", machine.getName().c_str());
    ImGui::Text("Type: Phase-Space | SAD: %.0f mm | Energy: %.0f MV",
                machine.getSAD(), geom.beamEnergyMV);
    ImGui::Text("PSF files: %d | MLC: %s (%d leaves)",
                geom.numPhaseSpaceFiles, geom.mlcType.c_str(), geom.numLeaves);

    // Show gantry angles from plan
    const auto& stfProps = m_state.plan->getStfProperties();
    ImGui::Text("Beams: %zu gantry angles", stfProps.gantryAngles.size());

    ImGui::Separator();

    renderLoadControls();

    if (m_state.phaseSpaceLoaded()) {
        ImGui::Separator();
        renderBeamVisibility();
        ImGui::Separator();
        renderStatistics();
        ImGui::Separator();
        renderVisualizationControls();
        ImGui::Separator();
        renderEnergyHistogram();
    }

    ImGui::End();
}

void PhaseSpacePanel::renderLoadControls() {
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Beam Configuration", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();

        // Show gantry angles from plan (read-only)
        const auto& stfProps = m_state.plan->getStfProperties();
        if (!stfProps.gantryAngles.empty()) {
            ImGui::Text("Gantry angles (from plan):");
            ImGui::SameLine();
            // Show first few angles
            std::string anglesStr;
            for (size_t i = 0; i < std::min(size_t(6), stfProps.gantryAngles.size()); ++i) {
                if (i > 0) anglesStr += ", ";
                anglesStr += std::to_string(static_cast<int>(stfProps.gantryAngles[i]));
            }
            if (stfProps.gantryAngles.size() > 6) {
                anglesStr += " ... (" + std::to_string(stfProps.gantryAngles.size()) + " total)";
            }
            ImGui::TextWrapped("%s", anglesStr.c_str());
        }

        ImGui::DragFloat("Collimator (deg)", &m_collimatorAngle, 1.0f, -180.0f, 180.0f, "%.1f");
        ImGui::DragFloat("Couch (deg)", &m_couchAngle, 1.0f, -90.0f, 90.0f, "%.1f");
        ImGui::Spacing();
        ImGui::DragInt("Max particles/beam (K)", &m_maxParticlesK, 100, 100, 50000);
        ImGui::DragInt("Viz sample/beam (K)", &m_vizSampleSizeK, 10, 10, 500);
        ImGui::Unindent();
    }

    ImGui::Spacing();

    if (m_isLoading) {
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Loading phase-space beams...");
        ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));

        if (m_loadDone) {
            m_loadThread.join();
            m_isLoading = false;
            m_loadDone = false;
            Logger::info("Phase-space multi-beam loading complete");
        }
    } else {
        const auto& stfProps = m_state.plan->getStfProperties();
        char buttonLabel[64];
        snprintf(buttonLabel, sizeof(buttonLabel), "Load Phase-Space (%zu beams)",
                 stfProps.gantryAngles.size());

        if (ImGui::Button(buttonLabel, ImVec2(-1, 30))) {
            m_state.resetPhaseSpace();
            m_isLoading = true;
            m_loadDone = false;
            m_selectedBeam = -1;

            m_loadThread = std::thread([this]() {
                auto start = std::chrono::steady_clock::now();

                try {
                    const auto& stfP = m_state.plan->getStfProperties();
                    const auto& gantryAngles = stfP.gantryAngles;

                    // Get isocenter from plan
                    std::array<double, 3> iso = {0.0, 0.0, 0.0};
                    if (!stfP.isoCenters.empty()) {
                        iso = stfP.isoCenters[0];
                    }

                    const int numBeams = static_cast<int>(gantryAngles.size());
                    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> sources(numBeams);

                    // Build all beams in parallel (OpenMP)
                    #pragma omp parallel for schedule(dynamic)
                    for (int i = 0; i < numBeams; ++i) {
                        auto source = std::make_shared<PhaseSpaceBeamSource>();
                        source->configure(m_state.plan->getMachine(),
                                           gantryAngles[i],
                                           m_collimatorAngle,
                                           m_couchAngle,
                                           iso);
                        source->build(static_cast<int64_t>(m_maxParticlesK) * 1000,
                                       static_cast<int64_t>(m_vizSampleSizeK) * 1000);
                        sources[i] = std::move(source);

                        Logger::debug("Phase-space beam " + std::to_string(i) +
                                      " (gantry=" + std::to_string(gantryAngles[i]) +
                                      ") loaded");
                    }

                    m_state.phaseSpaceSources = std::move(sources);

                    auto end = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration<double>(end - start).count();
                    m_statusMessage = "Loaded " + std::to_string(gantryAngles.size()) +
                                     " beams in " + std::to_string(elapsed).substr(0, 5) + "s";
                } catch (const std::exception& e) {
                    m_statusMessage = "Error: " + std::string(e.what());
                    Logger::error("Phase-space load failed: " + std::string(e.what()));
                }

                m_loadDone = true;
            });
        }
    }

    if (!m_statusMessage.empty()) {
        bool isError = m_statusMessage.find("Error") != std::string::npos;
        ImVec4 color = isError ? ImVec4(1.0f, 0.3f, 0.3f, 1.0f) : ImVec4(0.4f, 1.0f, 0.4f, 1.0f);
        ImGui::TextColored(color, "%s", m_statusMessage.c_str());
    }
}

void PhaseSpacePanel::renderBeamVisibility() {
    if (!m_phaseSpaceRenderer || m_phaseSpaceRenderer->getBeamCount() == 0) return;

    if (ImGui::CollapsingHeader("Beam Visibility", ImGuiTreeNodeFlags_DefaultOpen)) {
        // Show All / Hide All buttons
        if (ImGui::Button("Show All")) {
            m_phaseSpaceRenderer->setAllBeamsVisible(true);
        }
        ImGui::SameLine();
        if (ImGui::Button("Hide All")) {
            m_phaseSpaceRenderer->setAllBeamsVisible(false);
        }

        ImGui::Separator();

        // Per-beam checkboxes in a scrollable region
        ImGui::BeginChild("PhspBeamVisScroll", ImVec2(0, 200), true);
        for (size_t i = 0; i < m_phaseSpaceRenderer->getBeamCount(); ++i) {
            double gantryAngle = 0.0;
            if (i < m_state.phaseSpaceSources.size()) {
                gantryAngle = m_state.phaseSpaceSources[i]->getGantryAngle();
            }
            size_t numParticles = 0;
            if (i < m_state.phaseSpaceSources.size()) {
                numParticles = m_state.phaseSpaceSources[i]->getVisualizationSample().size();
            }

            bool visible = m_phaseSpaceRenderer->isBeamVisible(i);
            char label[64];
            snprintf(label, sizeof(label), "Beam %zu (%.0f deg, %zu pts)",
                     i, gantryAngle, numParticles);
            if (ImGui::Checkbox(label, &visible)) {
                m_phaseSpaceRenderer->setBeamVisible(i, visible);
            }
        }
        ImGui::EndChild();
    }
}

void PhaseSpacePanel::renderStatistics() {
    if (!m_state.phaseSpaceLoaded()) return;

    if (ImGui::CollapsingHeader("Verification Metrics", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();

        // Selector: aggregate or per-beam
        size_t nBeams = m_state.phaseSpaceSources.size();
        static std::vector<std::string> beamLabelStrs;
        static std::vector<const char*> items;
        beamLabelStrs.resize(nBeams);
        items.resize(nBeams + 1);
        items[0] = "All beams (aggregate)";
        for (size_t i = 0; i < nBeams; ++i) {
            char buf[48];
            snprintf(buf, sizeof(buf), "Beam %zu (%.0f deg)",
                     i, m_state.phaseSpaceSources[i]->getGantryAngle());
            beamLabelStrs[i] = buf;
            items[i + 1] = beamLabelStrs[i].c_str();
        }
        int comboIdx = m_selectedBeam + 1; // -1 → 0, 0 → 1, etc.
        if (ImGui::Combo("Stats for", &comboIdx, items.data(), static_cast<int>(nBeams + 1))) {
            m_selectedBeam = comboIdx - 1;
        }

        PhaseSpaceMetrics metrics;
        if (m_selectedBeam < 0) {
            // Aggregate across all beams
            int64_t totalCount = 0, photons = 0, electrons = 0, positrons = 0;
            double sumE = 0.0;
            double minE = 1e30, maxE = -1e30;
            for (const auto& src : m_state.phaseSpaceSources) {
                const auto& m = src->getMetrics();
                totalCount += m.totalCount;
                photons += m.photonCount;
                electrons += m.electronCount;
                positrons += m.positronCount;
                sumE += m.meanEnergy * m.totalCount;
                minE = std::min(minE, m.minEnergy);
                maxE = std::max(maxE, m.maxEnergy);
            }
            metrics.totalCount = totalCount;
            metrics.photonCount = photons;
            metrics.electronCount = electrons;
            metrics.positronCount = positrons;
            metrics.meanEnergy = totalCount > 0 ? sumE / totalCount : 0.0;
            metrics.minEnergy = minE;
            metrics.maxEnergy = maxE;
            // Take angular spread from first beam as representative
            if (!m_state.phaseSpaceSources.empty()) {
                metrics.angularSpreadU = m_state.phaseSpaceSources[0]->getMetrics().angularSpreadU;
                metrics.angularSpreadV = m_state.phaseSpaceSources[0]->getMetrics().angularSpreadV;
                metrics.xRange = m_state.phaseSpaceSources[0]->getMetrics().xRange;
                metrics.yRange = m_state.phaseSpaceSources[0]->getMetrics().yRange;
                metrics.zRange = m_state.phaseSpaceSources[0]->getMetrics().zRange;
            }
        } else if (static_cast<size_t>(m_selectedBeam) < m_state.phaseSpaceSources.size()) {
            metrics = m_state.phaseSpaceSources[m_selectedBeam]->getMetrics();
        }

        ImGui::Text("Total particles: %lld", static_cast<long long>(metrics.totalCount));
        ImGui::Text("  Photons:   %lld", static_cast<long long>(metrics.photonCount));
        ImGui::Text("  Electrons: %lld", static_cast<long long>(metrics.electronCount));
        ImGui::Text("  Positrons: %lld", static_cast<long long>(metrics.positronCount));
        ImGui::Spacing();
        ImGui::Text("Energy (MeV):");
        ImGui::Text("  Mean: %.4f", metrics.meanEnergy);
        ImGui::Text("  Range: [%.4f, %.4f]", metrics.minEnergy, metrics.maxEnergy);
        ImGui::Spacing();
        ImGui::Text("Angular spread:");
        ImGui::Text("  sigma_u: %.4f", metrics.angularSpreadU);
        ImGui::Text("  sigma_v: %.4f", metrics.angularSpreadV);
        ImGui::Spacing();
        ImGui::Text("Spatial extent (mm):");
        ImGui::Text("  X: [%.1f, %.1f]", metrics.xRange[0], metrics.xRange[1]);
        ImGui::Text("  Y: [%.1f, %.1f]", metrics.yRange[0], metrics.yRange[1]);
        ImGui::Text("  Z: [%.1f, %.1f]", metrics.zRange[0], metrics.zRange[1]);
        ImGui::Spacing();

        // Total viz sample count
        size_t totalViz = 0;
        for (const auto& src : m_state.phaseSpaceSources) {
            totalViz += src->getVisualizationSample().size();
        }
        ImGui::Text("Viz sample: %zu particles (%zu beams)",
                     totalViz, m_state.phaseSpaceSources.size());
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::Unindent();
    }
}

void PhaseSpacePanel::renderVisualizationControls() {
    if (!m_phaseSpaceRenderer) return;

    if (ImGui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();

        // Toggle particles
        bool showParts = m_phaseSpaceRenderer->isShowingParticles();
        if (ImGui::Checkbox("Show particles", &showParts)) {
            m_phaseSpaceRenderer->setShowParticles(showParts);
        }

        // Toggle direction vectors
        bool showDirs = m_phaseSpaceRenderer->isShowingDirections();
        if (ImGui::Checkbox("Show direction vectors", &showDirs)) {
            m_phaseSpaceRenderer->setShowDirections(showDirs);
        }

        // Point size
        float pointSize = m_phaseSpaceRenderer->getPointSize();
        if (ImGui::SliderFloat("Point size", &pointSize, 1.0f, 10.0f)) {
            m_phaseSpaceRenderer->setPointSize(pointSize);
        }

        // Direction scale
        if (showDirs) {
            float dirScale = m_phaseSpaceRenderer->getDirectionScale();
            if (ImGui::SliderFloat("Direction scale (mm)", &dirScale, 1.0f, 1000.0f)) {
                m_phaseSpaceRenderer->setDirectionScale(dirScale);
            }
        }

        // Color mode
        static const char* colorModes[] = { "Energy gradient", "Particle type" };
        int currentMode = static_cast<int>(m_phaseSpaceRenderer->getColorMode());
        if (ImGui::Combo("Color mode", &currentMode, colorModes, IM_ARRAYSIZE(colorModes))) {
            m_phaseSpaceRenderer->setColorMode(static_cast<PhaseSpaceColorMode>(currentMode));
        }

        ImGui::Unindent();
    }
}

void PhaseSpacePanel::renderEnergyHistogram() {
    if (!m_state.phaseSpaceLoaded()) return;

    if (ImGui::CollapsingHeader("Energy Distribution")) {
        // Use selected beam or first beam for histogram
        int beamIdx = m_selectedBeam >= 0 ? m_selectedBeam : 0;
        if (static_cast<size_t>(beamIdx) >= m_state.phaseSpaceSources.size()) return;

        auto histogram = m_state.phaseSpaceSources[beamIdx]->computeEnergyHistogram(40);
        if (!histogram.empty()) {
            // Find max count for normalization
            int64_t maxCount = 0;
            for (const auto& [energy, count] : histogram) {
                maxCount = std::max(maxCount, count);
            }

            if (maxCount > 0) {
                // Convert to float array for ImGui::PlotHistogram
                std::vector<float> values;
                values.reserve(histogram.size());
                for (const auto& [energy, count] : histogram) {
                    values.push_back(static_cast<float>(count));
                }

                char overlay[64];
                snprintf(overlay, sizeof(overlay), "Beam %d: %.3f - %.3f MeV",
                         beamIdx, histogram.front().first, histogram.back().first);

                ImGui::PlotHistogram("##EnergyHist", values.data(),
                                     static_cast<int>(values.size()),
                                     0, overlay, 0.0f,
                                     static_cast<float>(maxCount),
                                     ImVec2(-1, 120));
            }
        }
    }
}

} // namespace optirad
