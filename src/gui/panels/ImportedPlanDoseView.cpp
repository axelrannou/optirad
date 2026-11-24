#include "ImportedPlanDoseView.hpp"
#include "../Layout.hpp"
#include "Theme.hpp"
#include "core/workflow/ImportedPlanConverter.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <chrono>
#include <cstdio>

namespace optirad {

ImportedPlanDoseView::ImportedPlanDoseView(GuiAppState& state) : m_state(state) {}

ImportedPlanDoseView::~ImportedPlanDoseView() {
    if (m_thread.joinable()) m_thread.join();
}

void ImportedPlanDoseView::installImported(GuiAppState& state, int importedDoseId) {
    const ImportedPlan* imported = state.patientData ? state.patientData->getImportedPlan() : nullptr;
    if (!imported || imported->beams.empty()) return;
    try {
        installPrepared(state, importedDoseId, prepareImported(state.patientData, *imported));
    } catch (const std::exception& e) {
        Logger::error("RTPLAN plan could not be built: " + std::string(e.what()));
    }
    state.syncSelectedDose();
}

ImportedPlanDoseView::PreparedPlan ImportedPlanDoseView::prepareImported(std::shared_ptr<PatientData> patient,
                                                                         const ImportedPlan& imported) {
    const auto start = std::chrono::steady_clock::now();
    PreparedPlan prepared;
    PlanConfig config;
    ImportedPlanConverter::applyToPlanConfig(imported, config);
    if (const auto& rtDoseGrid = patient->getImportedDoseGridResolution())
        config.doseGridResolution = *rtDoseGrid;
    auto& built = prepared.built;
    built = PlanBuilder::build(config, patient);
    for (const auto& w : ImportedPlanConverter::validate(imported, built.plan->getMachine()))
        Logger::warn("RTPLAN: " + w);

    auto& delivery = prepared.delivery;
    auto& segments = prepared.segments;
    if (imported.hasArcs()) {
        const auto mismatches = ImportedPlanConverter::arcPlanMismatches(imported, *built.plan);
        if (!mismatches.empty())
            throw std::runtime_error(mismatches.front());
        delivery = std::make_shared<ArcPipelineResult>(ArcPipeline::prepareSequence(
            ImportedPlanConverter::toArcSequence(imported, *built.plan), *built.plan, *patient));
    } else {
        segments = ImportedPlanConverter::toLeafSequences(imported, *built.plan);
        const size_t beams = built.stf ? built.stf->getCount() : 0;
        if (segments.size() != beams) {
            Logger::warn("RTPLAN: " + std::to_string(segments.size()) + " static beam(s) but " +
                         std::to_string(beams) + " plan beam(s); delivery not shown");
            segments.clear();
        }
    }

    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    Logger::info("RTPLAN '" + imported.label + "' prepared in " + std::to_string(s).substr(0, 5) + " s");
    return prepared;
}

void ImportedPlanDoseView::installPrepared(GuiAppState& state, int importedDoseId, PreparedPlan prepared) {
    auto& built = prepared.built;
    auto& delivery = prepared.delivery;
    auto& segments = prepared.segments;
    state.switchToPlan(nullptr);
    state.plan = built.plan;
    state.stfProps = built.stfProps;
    state.stf = built.stf;
    state.workingPlan = built.plan;
    state.importedPlan = built.plan;
    state.importedDelivery = delivery;
    state.importedSequences = segments;
    if (const DoseEntry* entry = state.doseStore.getEntry(state.doseStore.getSelectedIdx());
        entry && entry->id == importedDoseId) {
        if (delivery) {
            auto withDose = std::make_shared<ArcPipelineResult>(*delivery);
            withDose->dose = entry->dose;
            state.arcCache[importedDoseId] = {withDose, -1};
        }
        if (!segments.empty())
            state.seqCache[importedDoseId].sequences = segments;
        state.dosePlans[importedDoseId] = built.plan;
    }
    state.syncSelectedDose();
}

void ImportedPlanDoseView::render(const GuiAppState::DaoDoseSettings& dose) {
    const ImportedPlan* imported = m_state.patientData ? m_state.patientData->getImportedPlan() : nullptr;
    if (!imported || imported->beams.empty()) return;

    if (!m_state.importedPlan) {
        ImGui::TextColored(getThemeColors().failText, "No plan could be built from the RTPLAN (see the log).");
    } else if (!imported->hasArcs()) {
        ImGui::TextDisabled("%s", m_state.importedSequences.empty()
            ? "Static-beam RTPLAN: its segments could not be matched to the plan beams (see the log)."
            : "Static-beam RTPLAN: segments shown in Leaf Sequencing, BEV and playback; dose recompute from "
              "segments is not supported yet.");
    } else {
        const auto& rtDoseGrid = m_state.patientData->getImportedDoseGridResolution();
        if (rtDoseGrid) {
            char label[128];
            snprintf(label, sizeof(label), "RTDOSE grid (%.2f x %.2f x %.2f mm)",
                     (*rtDoseGrid)[0], (*rtDoseGrid)[1], (*rtDoseGrid)[2]);
            ImGui::Checkbox(label, &m_useRtDoseGrid);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Off: the dose grid resolution of the DAO dose settings.");
        }

        const bool busy = m_running || m_state.taskRunning || !m_state.importedDelivery;
        const char* label = "Compute dose from RTPLAN leaves";
        if (busy) ImGui::BeginDisabled();
        if (ImGui::Button(label, layout::fullWidthButtonSize(label, 30))) start(dose);
        if (busy) ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
            ImGui::SetTooltip("New dose: the delivered dose of the RTPLAN's own leaf positions and MU,\n"
                              "without optimization, computed with the %s machine model.",
                              m_state.importedPlan->getMachine().getName().c_str());
    }

    if (m_running) {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            ImGui::TextColored(getThemeColors().progressText, "%s...", m_stage.c_str());
        }
        ImGui::ProgressBar(-1.0f * static_cast<float>(ImGui::GetTime()), ImVec2(-1, 0));
        if (m_done) {
            m_thread.join();
            m_running = false;
            m_done = false;
            m_state.taskRunning = false;
            if (m_pending && m_pendingPlan == m_state.importedPlan) {
                const std::string name = "RTPLAN " + (imported->label.empty() ? std::string("dose") : imported->label) +
                                         " (recomputed)";
                const int doseId = m_state.doseStore.addDose(name, m_pending->dose, m_pendingGrid);
                m_state.dosePlans[doseId] = m_pendingPlan;
                m_state.arcCache[doseId] = {m_pending, -1};
                m_state.syncSelectedDose();
                m_state.optimizationJustFinished = true;
            }
            m_pending.reset();
            m_pendingGrid.reset();
            m_pendingPlan.reset();
        }
    } else if (!m_status.empty()) {
        ImGui::TextColored(m_failed ? getThemeColors().failText : getThemeColors().passText, "%s", m_status.c_str());
    }
}

void ImportedPlanDoseView::start(const GuiAppState::DaoDoseSettings& dose) {
    if (m_thread.joinable()) m_thread.join();
    m_running = true;
    m_done = false;
    m_failed = false;
    m_state.taskRunning = true;
    m_status.clear();
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_stage = "Preparing";
    }

    ArcOptimizationConfig cfg;
    cfg.computeBackend = dose.backend;
    cfg.gpuMemoryBudgetMiB = dose.gpuMemoryBudgetMiB;
    cfg.engineOptions.numThreads = dose.numThreads;
    cfg.engineOptions.excludeExternal = dose.excludeExternal;
    cfg.fluencePixelMm = dose.fluencePixelMm;

    const auto& rtDoseGrid = m_state.patientData->getImportedDoseGridResolution();
    const std::array<double, 3> resolution = m_useRtDoseGrid && rtDoseGrid ? *rtDoseGrid : dose.resolution;
    auto patient = m_state.patientData;
    m_pendingPlan = m_state.importedPlan;
    auto delivery = m_state.importedDelivery;

    m_thread = std::thread([this, cfg, resolution, patient, delivery]() {
        try {
            const Plan& plan = *m_pendingPlan;
            m_pendingGrid = std::make_shared<Grid>(Grid::createDoseGrid(patient->getGrid(), resolution));
            auto result = std::make_shared<ArcPipelineResult>(*delivery);
            ArcPipeline::computeDeliveredDose(*result, plan, *patient, *m_pendingGrid, cfg, plan.getPrescribedDose(),
                                              [this](const std::string& stage) {
                                                  std::lock_guard<std::mutex> lock(m_mutex);
                                                  m_stage = stage;
                                              });
            const auto& seq = result->sequence;
            char buf[256];
            snprintf(buf, sizeof(buf), "%zu control points | %.0f MU/fraction | %.0f s delivery",
                     seq.controlPoints.size(), seq.totalMU * plan.machineMuPerWeight(), seq.totalTimeSec);
            m_status = buf;
            m_pending = std::move(result);
        } catch (const std::exception& e) {
            Logger::error("Dose from RTPLAN failed: " + std::string(e.what()));
            m_status = "Error: " + std::string(e.what());
            m_failed = true;
            m_pending.reset();
        }
        m_done = true;
    });
}

} // namespace optirad
