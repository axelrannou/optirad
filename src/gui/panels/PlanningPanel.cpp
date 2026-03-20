#include "PlanningPanel.hpp"
#include "Theme.hpp"
#include "io/MachineLoader.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "dose/DijSerializer.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>
#include <chrono>
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

        // --- Parse gantry angles ---
        if (m_gantryMode == 0) {
            stfProps.setGantryAngles(m_gantryStart, m_gantryStep, m_gantryStop);
        } else {
            std::vector<double> gList;
            { std::istringstream ss(m_gantryListBuf); double v; while (ss >> v) gList.push_back(v); }
            stfProps.setGantryAngles(gList);
        }

        // --- Parse couch angles ---
        if (m_couchMode == 0) {
            if (m_couchStep > 0.0f) {
                stfProps.setCouchAngles(
                    static_cast<double>(m_couchStart),
                    static_cast<double>(m_couchStep),
                    static_cast<double>(m_couchStop));
            } else {
                stfProps.setUniformCouchAngle(static_cast<double>(m_couchStart));
            }
        } else {
            std::vector<double> cList;
            { std::istringstream ss(m_couchListBuf); double v; while (ss >> v) cList.push_back(v); }
            stfProps.setCouchAngles(cList);
        }

        stfProps.ensureConsistentAngles();
        stfProps.bixelWidth = m_bixelWidth;

        // Compute isocenter from target structures
        auto iso = plan->computeIsoCenter();
        stfProps.setUniformIsoCenter(iso);

        plan->setStfProperties(stfProps);

        m_state.plan = plan;
        m_planCreated = true;

        Logger::info("Plan created: " + std::to_string(stfProps.numOfBeams) +
            " beams, bixelWidth=" + std::to_string(m_bixelWidth) + "mm");

        // Automatically generate STF for generic (non-phase-space) machines
        if (!kIsPhaseSpace[m_machineIdx]) {
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

                double bixelWidth = stfP.bixelWidth;

                PhotonIMRTStfGenerator generator(0.0, 360.0, 360.0, bixelWidth, iso);
                generator.setMachine(m_state.plan->getMachine());
                generator.setRadiationMode(radiationMode);
                generator.setGantryAngles(stfP.gantryAngles);
                generator.setCouchAngles(stfP.couchAngles);

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
                m_stfStatusMessage = "STF generated in " + std::to_string(elapsed).substr(0, 5) + "s";

                m_stfGenerationDone = true;
            });
        }
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
                Logger::info("Dij calculation complete");
            }
        } else {
            if (!canCalcDose) ImGui::BeginDisabled();

            if (ImGui::Button("Calculate Dij", ImVec2(-1, 30))) {
                m_state.resetDij();
                m_isCalculatingDij = true;
                m_dijCalcDone = false;
                m_state.cancelFlag = false;
                m_dijCurrentBeam = 0;
                m_dijTotalBeams = 0;

                m_dijThread = std::thread([this]() {
                    auto start = std::chrono::steady_clock::now();

                    try {
                        // Update plan's dose grid resolution
                        m_state.plan->setDoseGridResolution({
                            static_cast<double>(m_doseResolution[0]),
                            static_cast<double>(m_doseResolution[1]),
                            static_cast<double>(m_doseResolution[2])});

                        // Create dose grid
                        const auto& ctGrid = m_state.patientData->getGrid();
                        auto doseGrid = std::make_shared<Grid>(
                            Grid::createDoseGrid(ctGrid, {
                                static_cast<double>(m_doseResolution[0]),
                                static_cast<double>(m_doseResolution[1]),
                                static_cast<double>(m_doseResolution[2])}));
                        m_state.doseGrid = doseGrid;

                        auto dims = doseGrid->getDimensions();
                        Logger::info("Dose grid: " + std::to_string(dims[0]) + "x" +
                            std::to_string(dims[1]) + "x" + std::to_string(dims[2]));

                        // Check cache
                        std::string patientName = "unknown";
                        if (m_state.patientData->getPatient()) {
                            patientName = m_state.patientData->getPatient()->getName();
                        }
                        std::string cacheFile = DijSerializer::getCacheDir() + "/" +
                            DijSerializer::buildCacheKey(
                                patientName,
                                static_cast<int>(m_state.stf->getCount()),
                                m_state.plan->getStfProperties().bixelWidth,
                                static_cast<double>(m_doseResolution[0]));

                        if (DijSerializer::exists(cacheFile)) {
                            Logger::info("Loading Dij from cache: " + cacheFile);
                            auto loaded = DijSerializer::load(cacheFile);
                            m_state.dij = std::make_shared<DoseInfluenceMatrix>(std::move(loaded));
                            m_dijStatusMessage = "Loaded from cache";
                        } else {
                            // Create dose engine and compute
                            auto engine = DoseEngineFactory::create("PencilBeam");
                            engine->setCancelFlag(&m_state.cancelFlag);
                            engine->setProgressCallback(
                                [this](int current, int total, const std::string&) {
                                    m_dijCurrentBeam = current;
                                    m_dijTotalBeams = total;
                                });

                            // Apply user options (thresholds, threads)
                            DoseCalcOptions opts;
                            opts.absoluteThreshold = static_cast<double>(m_absoluteThreshold);
                            opts.relativeThreshold = static_cast<double>(m_relativeThreshold) / 100.0;
                            opts.numThreads = m_numThreads;
                            engine->setOptions(opts);

                            auto dij = engine->calculateDij(
                                *m_state.plan, *m_state.stf, *m_state.patientData, *doseGrid);
                            m_state.dij = std::make_shared<DoseInfluenceMatrix>(std::move(dij));

                            // Save to cache
                            DijSerializer::save(*m_state.dij, cacheFile);

                            auto end = std::chrono::steady_clock::now();
                            double elapsed = std::chrono::duration<double>(end - start).count();
                            m_dijStatusMessage = "Computed in " +
                                std::to_string(elapsed).substr(0, 5) + "s" +
                                " (nnz=" + std::to_string(m_state.dij->getNumNonZeros()) + ")";
                        }
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
