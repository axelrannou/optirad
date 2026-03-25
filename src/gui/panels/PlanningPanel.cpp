#include "PlanningPanel.hpp"
#include "Theme.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "core/workflow/DoseCalculationPipeline.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>
#include <sstream>

namespace optirad {

static const char* kRadiationModes[] = { "photons" };
static const char* kMachines[] = { "Generic", "Varian_TrueBeam6MV" };
static const bool kIsPhaseSpace[] = { false, true };

// Helper callback for ImGui InputText with std::string
static int InputTextCallback(ImGuiInputTextCallbackData* data) {
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        std::string* str = static_cast<std::string*>(data->UserData);
        str->resize(data->BufTextLen);
        data->Buf = &(*str)[0];
    }
    return 0;
}

// Helper wrapper for InputTextMultiline with std::string
static bool InputTextMultilineString(const char* label, std::string* str, const ImVec2& size = ImVec2(0, 0), ImGuiInputTextFlags flags = 0) {
    flags |= ImGuiInputTextFlags_CallbackResize;
    str->reserve(256);  // Initial capacity
    return ImGui::InputTextMultiline(label, &(*str)[0], str->capacity() + 1, size, flags, InputTextCallback, str);
}

PlanningPanel::PlanningPanel(GuiAppState& state) : m_state(state) {}

PlanningPanel::~PlanningPanel() {
    if (m_stfThread.joinable()) m_stfThread.join();
    if (m_dijThread.joinable()) m_dijThread.join();
}

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
        const auto& col = getThemeColors();
        if (kIsPhaseSpace[m_machineIdx]) {
            ImGui::TextColored(col.phaseSpace, "Phase-space (IAEA PSF)");
        } else {
            ImGui::TextColored(col.generic, "Generic (pencil-beam)");
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
        ImGui::RadioButton("Range##gantry", &m_gantryMode, 0); ImGui::SameLine();
        ImGui::RadioButton("List##gantry", &m_gantryMode, 1);

        if (m_gantryMode == 0) {
            ImGui::InputFloat("Start (deg)##gantry", &m_gantryStart, 0.0f, 0.0f, "%.1f");
            ImGui::InputFloat("Step (deg)##gantry", &m_gantryStep, 0.0f, 0.0f, "%.1f");
            ImGui::InputFloat("Stop (deg)##gantry", &m_gantryStop, 0.0f, 0.0f, "%.1f");
            int numGantry = 0;
            if (m_gantryStep > 0.0f && m_gantryStop > m_gantryStart) {
                numGantry = static_cast<int>(std::ceil((m_gantryStop - m_gantryStart) / m_gantryStep));
            }
            ImGui::Text("Gantry entries: %d", numGantry);
        } else {
            InputTextMultilineString("##gantryList", &m_gantryListBuf, ImVec2(-1, ImGui::GetTextLineHeight() * 4));
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Space-separated gantry angles in degrees.\nExample: 0 90 180 270");
            // Count entries
            int cnt = 0;
            { std::istringstream ss(m_gantryListBuf); double v; while (ss >> v) ++cnt; }
            ImGui::Text("Gantry entries: %d", cnt);
        }
        ImGui::Unindent();
    }

    // ── Couch Angles ──
    if (ImGui::CollapsingHeader("Couch Angles", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        ImGui::RadioButton("Range##couch", &m_couchMode, 0); ImGui::SameLine();
        ImGui::RadioButton("List##couch", &m_couchMode, 1);

        if (m_couchMode == 0) {
            ImGui::InputFloat("Start (deg)##couch", &m_couchStart, 0.0f, 0.0f, "%.1f");
            ImGui::InputFloat("Step (deg)##couch", &m_couchStep, 0.0f, 0.0f, "%.1f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Set to 0 for a single couch angle (= Start) on every beam.\n"
                                  "Step > 0 creates a multi-arc plan: each couch angle\n"
                                  "sweeps through ALL gantry angles (Cartesian product).\n"
                                  "Total beams = gantry_entries x couch_entries.");
            ImGui::InputFloat("Stop (deg)##couch", &m_couchStop, 0.0f, 0.0f, "%.1f");

            if (m_couchStep <= 0.0f) {
                ImGui::Text("Single couch angle: %.1f deg (all beams)", m_couchStart);
            } else {
                int numCouch = 0;
                if (m_couchStop > m_couchStart) {
                    numCouch = static_cast<int>(std::ceil((m_couchStop - m_couchStart) / m_couchStep));
                }
                ImGui::Text("Couch entries: %d", numCouch);
            }
        } else {
            InputTextMultilineString("##couchList", &m_couchListBuf, ImVec2(-1, ImGui::GetTextLineHeight() * 4));
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Space-separated couch angles in degrees (paired 1:1 with gantry).\nExample: 0 10 20 30");
            int cnt = 0;
            { std::istringstream ss(m_couchListBuf); double v; while (ss >> v) ++cnt; }
            ImGui::Text("Couch entries: %d", cnt);
        }
        ImGui::Unindent();
    }

    // ── Beam count ──
    {
        int numGantry = 0;
        if (m_gantryMode == 0) {
            if (m_gantryStep > 0.0f && m_gantryStop > m_gantryStart)
                numGantry = static_cast<int>(std::ceil((m_gantryStop - m_gantryStart) / m_gantryStep));
        } else {
            std::istringstream ss(m_gantryListBuf); double v; while (ss >> v) ++numGantry;
        }

        int numBeams = numGantry;
        if (m_couchMode == 0 && m_couchStep > 0.0f && m_couchStop > m_couchStart) {
            int numCouch = static_cast<int>(std::ceil((m_couchStop - m_couchStart) / m_couchStep));
            if (numCouch > 1) {
                numBeams = numGantry * numCouch;
                ImGui::Text("Number of beams: %d  (%d gantry x %d arcs)",
                            numBeams, numGantry, numCouch);
            } else {
                ImGui::Text("Number of beams: %d", numBeams);
            }
        } else if (m_couchMode == 1) {
            // Paired 1:1 — beam count = max of gantry/couch lists
            int numCouchList = 0;
            { std::istringstream ss(m_couchListBuf); double v; while (ss >> v) ++numCouchList; }
            numBeams = std::max(numGantry, numCouchList);
            ImGui::Text("Number of beams: %d  (paired 1:1)", numBeams);
        } else {
            ImGui::Text("Number of beams: %d", numBeams);
        }
    }

    // ── Bixel Width (only for generic machines) ──
    if (!kIsPhaseSpace[m_machineIdx]) {
        if (ImGui::CollapsingHeader("Bixel Width", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Indent();
            ImGui::InputFloat("Width (mm)", &m_bixelWidth, 0.0f, 0.0f, "%.1f");
            ImGui::Unindent();
        }
    }

    ImGui::Separator();
    ImGui::Spacing();

    // ── Validation ──
    bool paramsValid = false;
    if (m_gantryMode == 0) {
        int numGantryV = 0;
        if (m_gantryStep > 0.0f && m_gantryStop > m_gantryStart)
            numGantryV = static_cast<int>(std::ceil((m_gantryStop - m_gantryStart) / m_gantryStep));
        paramsValid = (numGantryV > 0);
    } else {
        int cnt = 0;
        { std::istringstream ss(m_gantryListBuf); double v; while (ss >> v) ++cnt; }
        paramsValid = (cnt > 0);
    }
    paramsValid = paramsValid &&
                  (kIsPhaseSpace[m_machineIdx] || m_bixelWidth > 0.0f) &&
                  (m_numFractions > 0);

    if (!paramsValid) {
        ImGui::TextColored(getThemeColors().failText, "Fix parameter errors above.");
    }

    // ── Create Plan button ──
    if (!paramsValid) ImGui::BeginDisabled();

    if (ImGui::Button("Create Plan", ImVec2(-1, 30))) {
        // Reset downstream state
        m_state.resetPlan();

        // Build config from UI values
        PlanConfig config;
        config.radiationMode = kRadiationModes[m_radiationModeIdx];
        config.machineName = kMachines[m_machineIdx];
        config.numFractions = m_numFractions;
        config.bixelWidth = m_bixelWidth;

        if (m_gantryMode == 0) {
            config.gantryStart = m_gantryStart;
            config.gantryStep = m_gantryStep;
            config.gantryStop = m_gantryStop;
        } else {
            std::istringstream ss(m_gantryListBuf);
            double v;
            while (ss >> v) config.gantryAngles.push_back(v);
        }

        if (m_couchMode == 0) {
            config.couchStart = m_couchStart;
            config.couchStep = m_couchStep;
            config.couchStop = m_couchStop;
        } else {
            std::istringstream ss(m_couchListBuf);
            double v;
            while (ss >> v) config.couchAngles.push_back(v);
        }

        m_isGeneratingStf = true;
        m_stfGenerationDone = false;
        m_state.taskRunning = true;

        m_stfThread = std::thread([this, config]() {
            auto start = std::chrono::steady_clock::now();
            try {
                auto result = PlanBuilder::build(config, m_state.patientData);
                m_state.plan = result.plan;
                if (result.stfProps) m_state.stfProps = result.stfProps;
                if (result.stf) m_state.stf = result.stf;
                m_planCreated = true;

                const auto& stfP = m_state.plan->getStfProperties();
                Logger::info("Plan created: " + std::to_string(stfP.numOfBeams) +
                    " beams, bixelWidth=" + std::to_string(config.bixelWidth) + "mm");

                auto end = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(end - start).count();
                m_stfStatusMessage = "Plan created in " + std::to_string(elapsed).substr(0, 5) + "s";
            } catch (const std::exception& e) {
                Logger::error("Plan creation failed: " + std::string(e.what()));
                m_stfStatusMessage = "Error: " + std::string(e.what());
            }
            m_stfGenerationDone = true;
        });
    }

    if (!paramsValid) ImGui::EndDisabled();

    // Handle async STF generation completion
    if (m_isGeneratingStf) {
        ImGui::TextColored(getThemeColors().progressText, "Creating plan...");
        ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));

        if (m_stfGenerationDone) {
            m_stfThread.join();
            m_isGeneratingStf = false;
            m_stfGenerationDone = false;
            m_state.taskRunning = false;
        }
    }

    // Show plan status
    if (m_state.planCreated()) {
        ImGui::Spacing();
        if (m_state.isPhaseSpaceMachine()) {
            ImGui::TextColored(getThemeColors().passText, "Plan created.");
            ImGui::TextColored(getThemeColors().phaseSpace,
                "Phase-space machine: use the Phase Space panel to load beam data.");
        } else if (m_state.stfGenerated()) {
            ImGui::TextColored(getThemeColors().passText, "Plan created.");
            if (!m_stfStatusMessage.empty()) {
                ImGui::TextColored(getThemeColors().passText, "%s", m_stfStatusMessage.c_str());
            }
            ImGui::Text("Beams: %zu | Rays: %zu | Bixels: %zu",
                        m_state.stf->getCount(),
                        m_state.stf->getTotalNumOfRays(),
                        m_state.stf->getTotalNumOfBixels());
        } //else {
        //     ImGui::TextColored(getThemeColors().passText, "Plan created.");
        // }

        const auto& stfProps = m_state.plan->getStfProperties();
        ImGui::Text("  Beams: %zu", stfProps.gantryAngles.size());
        ImGui::Text("  Iso: (%.1f, %.1f, %.1f) mm",
                    stfProps.isoCenters[0][0],
                    stfProps.isoCenters[0][1],
                    stfProps.isoCenters[0][2]);
    }

    // ── Dose Calculation ──
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Dose Calculation", ImGuiTreeNodeFlags_DefaultOpen)) {
        bool canCalcDose = m_state.stfGenerated() && !m_state.isPhaseSpaceMachine();

        if (!m_state.stfGenerated()) {
            ImGui::TextDisabled("Create plan first to calculate dose.");
        }

        // Dose grid resolution
        if (canCalcDose && !m_isCalculatingDij) {
            ImGui::Text("Dose Grid Resolution:");
            ImGui::SetNextItemWidth(95.0f);
            ImGui::InputFloat("##doseRes0", &m_doseResolution[0], 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(95.0f);
            ImGui::InputFloat("##doseRes1", &m_doseResolution[1], 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(95.0f);
            ImGui::InputFloat("##doseRes2", &m_doseResolution[2], 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::TextUnformatted("mm");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Text("Memory / Speed Options:");

            ImGui::SetNextItemWidth(120);
            ImGui::InputFloat("Relative Threshold (%%)", &m_relativeThreshold, 0.0f, 0.0f, "%.1f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Discard entries below this %% of each bixel's max dose.\n"
                                  "Higher = less RAM, slightly less accurate.\n"
                                  "Recommended: 1%%. Set 0 to disable.");

            ImGui::SetNextItemWidth(120);
            ImGui::InputFloat("Absolute Threshold", &m_absoluteThreshold, 0.0f, 0.0f, "%.1f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Discard dose entries below this value (Gy).\n"
                                  "Removes noise from far-field kernel tails.\n"
                                  "Recommended: 1e-6.");
            if (m_absoluteThreshold < 0.0f) m_absoluteThreshold = 0.0f;

            ImGui::SetNextItemWidth(120);
            ImGui::InputInt("Threads (0=all)", &m_numThreads);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Number of OpenMP threads.\n"
                                  "0 = use all available CPU cores.");
            if (m_numThreads < 0) m_numThreads = 0;
        }

        if (m_isCalculatingDij) {
            ImGui::TextColored(getThemeColors().progressText, "Calculating Dij...");
            if (m_dijTotalBeams > 0) {
                float progress = static_cast<float>(m_dijCurrentBeam) / static_cast<float>(m_dijTotalBeams);
                ImGui::ProgressBar(progress, ImVec2(-1, 0),
                    (std::to_string(m_dijCurrentBeam) + "/" + std::to_string(m_dijTotalBeams) + " beams").c_str());
            } else {
                ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));
            }

            if (ImGui::Button("Cancel", ImVec2(-1, 24))) {
                m_state.cancelFlag = true;
            }

            if (m_dijCalcDone) {
                m_dijThread.join();
                m_isCalculatingDij = false;
                m_dijCalcDone = false;
                m_state.cancelFlag = false;
                m_state.taskRunning = false;

                Logger::info("Dij calculation complete");
            }
        } else {
            if (!canCalcDose) ImGui::BeginDisabled();

            if (ImGui::Button("Calculate Dij", ImVec2(-1, 30))) {
                m_state.resetDij();
                m_isCalculatingDij = true;
                m_dijCalcDone = false;
                m_state.cancelFlag = false;
                m_state.taskRunning = true;
                m_dijCurrentBeam = 0;
                m_dijTotalBeams = 0;

                m_dijThread = std::thread([this]() {
                    try {
                        m_state.plan->setDoseGridResolution({
                            static_cast<double>(m_doseResolution[0]),
                            static_cast<double>(m_doseResolution[1]),
                            static_cast<double>(m_doseResolution[2])});

                        DoseCalcPipelineOptions opts;
                        opts.resolution = {
                            static_cast<double>(m_doseResolution[0]),
                            static_cast<double>(m_doseResolution[1]),
                            static_cast<double>(m_doseResolution[2])};
                        opts.absoluteThreshold = static_cast<double>(m_absoluteThreshold);
                        opts.relativeThreshold = static_cast<double>(m_relativeThreshold) / 100.0;
                        opts.numThreads = m_numThreads;

                        auto result = DoseCalculationPipeline::run(
                            *m_state.plan, *m_state.stf, *m_state.patientData, opts,
                            [this](int current, int total) {
                                m_dijCurrentBeam = current;
                                m_dijTotalBeams = total;
                            },
                            &m_state.cancelFlag);

                        m_state.computeGrid = result.doseGrid;
                        m_state.dij = result.dij;
                        m_dijStatusMessage = result.cacheHit ? "Loaded from cache" :
                            "Computed (nnz=" + std::to_string(result.dij->getNumNonZeros()) + ")";
                    } catch (const std::exception& e) {
                        Logger::error("Dij calculation failed: " + std::string(e.what()));
                        m_dijStatusMessage = "Error: " + std::string(e.what());
                    }

                    m_dijCalcDone = true;
                });
            }

            if (!canCalcDose) ImGui::EndDisabled();
        }

        // Show Dij status
        if (m_state.dijComputed()) {
            if (!m_dijStatusMessage.empty()) {
                ImGui::TextColored(getThemeColors().passText, "%s", m_dijStatusMessage.c_str());
            }
            ImGui::Text("Dij: %zu voxels x %zu bixels (nnz: %zu)",
                        m_state.dij->getNumVoxels(),
                        m_state.dij->getNumBixels(),
                        m_state.dij->getNumNonZeros());
        }
    }

    ImGui::End();
}

} // namespace optirad
