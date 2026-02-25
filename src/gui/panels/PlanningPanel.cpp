#include "PlanningPanel.hpp"
#include "io/MachineLoader.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "dose/DijSerializer.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <cmath>
#include <chrono>

namespace optirad {

static const char* kRadiationModes[] = { "photons" };
static const char* kMachines[] = { "Generic", "Varian_TrueBeam6MV" };
static const bool kIsPhaseSpace[] = { false, true };

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

    // ── Dose Calculation ──
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Dose Calculation", ImGuiTreeNodeFlags_DefaultOpen)) {
        bool canCalcDose = m_state.stfGenerated() && !m_state.isPhaseSpaceMachine();

        if (!m_state.stfGenerated()) {
            ImGui::TextDisabled("Generate STF first to calculate dose.");
        }

        // Dose grid resolution
        if (canCalcDose && !m_isCalculatingDij) {
            ImGui::Text("Dose Grid Resolution (mm):");
            ImGui::DragFloat3("##doseRes", m_doseResolution, 0.1f, 0.5f, 10.0f, "%.1f");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Text("Memory / Speed Options:");

            ImGui::SetNextItemWidth(120);
            ImGui::DragFloat("Relative Threshold (%%)", &m_relativeThreshold, 0.1f, 0.0f, 50.0f, "%.1f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Discard entries below this %% of each bixel's max dose.\n"
                                  "Higher = less RAM, slightly less accurate.\n"
                                  "Recommended: 1%%. Set 0 to disable.");

            ImGui::SetNextItemWidth(120);
            ImGui::InputFloat("Absolute Threshold", &m_absoluteThreshold, 0.0f, 0.0f, "%.1e");
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
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Calculating Dij...");
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
                ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "%s", m_dijStatusMessage.c_str());
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
