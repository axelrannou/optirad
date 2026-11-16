#include "DaoOptimizationView.hpp"
#include "Theme.hpp"
#include "optimization/ObjectiveBuilder.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

namespace optirad {

DaoOptimizationView::DaoOptimizationView(GuiAppState& state) : m_state(state) {}

DaoOptimizationView::~DaoOptimizationView() {
    if (m_thread.joinable()) m_thread.join();
}

void DaoOptimizationView::renderSettings() {
    ImGui::InputInt("Iterations per resolution", &m_iterations);
    m_iterations = std::clamp(m_iterations, 10, 5000);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("The last (finest) resolution gets twice as many.");
    ImGui::InputFloat("First CP spacing (deg)", &m_coarsestSpacing, 0.0f, 0.0f, "%.1f");
    m_coarsestSpacing = std::max(m_coarsestSpacing, 1.0f);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Control point spacing of the first resolution; halved until the plan's spacing\n"
                          "(progressive resolution, as in Eclipse).");
    ImGui::InputFloat("MU penalty", &m_muPenalty, 0.05f, 0.25f, "%.2f");
    m_muPenalty = std::max(m_muPenalty, 0.0f);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Cost of the total MU as a fraction of the objective at the second resolution. 0 = none.");
    ImGui::InputFloat("MU smoothing", &m_muSmoothing, 0.01f, 0.1f, "%.3f");
    m_muSmoothing = std::max(m_muSmoothing, 0.0f);
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Penalty on MU differences between neighbouring control points.");
    ImGui::Checkbox("Integrate leaf motion", &m_integrateMotion);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Dose of a control point = average over the leaf motion in its slice.");
}

void DaoOptimizationView::start(std::vector<ObjectiveSpec> specs, const OptimizationConfig& config) {
    if (m_thread.joinable()) m_thread.join();
    const auto& dose = m_state.daoDose;
    m_state.plan->setPrescribedDose(config.targetDose);
    m_state.plan->setDoseGridResolution(dose.resolution);
    m_state.computeGrid = std::make_shared<Grid>(Grid::createDoseGrid(m_state.patientData->getGrid(), dose.resolution));
    m_state.computeBackend = dose.backend;
    m_state.gpuMemoryBudgetMiB = dose.gpuMemoryBudgetMiB;
    m_state.optimizationSpecs = specs;
    m_state.optimizationConfig = config;
    m_state.resetLeafSequence();

    ArcOptimizationConfig cfg;
    cfg.base = config;
    cfg.iterationsPerResolution = m_iterations;
    cfg.coarsestSpacingDeg = static_cast<double>(m_coarsestSpacing);
    cfg.muPenalty = static_cast<double>(m_muPenalty);
    cfg.muSmoothingWeight = static_cast<double>(m_muSmoothing);
    cfg.integrateLeafMotion = m_integrateMotion;
    const double resolution = m_state.plan->getMachine().getGeometry().leafPositionResolution;
    if (resolution > 0.0) cfg.seed.leafPositionResolution = resolution;
    cfg.computeBackend = dose.backend;
    cfg.gpuMemoryBudgetMiB = dose.gpuMemoryBudgetMiB;
    cfg.engineOptions.numThreads = dose.numThreads;
    cfg.engineOptions.excludeExternal = dose.excludeExternal;
    cfg.fluencePixelMm = dose.fluencePixelMm;

    m_running = true;
    m_done = false;
    m_state.taskRunning = true;
    m_status.clear();
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_progress.clear();
        m_stage = "Preparing DAO";
    }
    m_displayHistory.clear();
    m_progressStage = 0;

    m_thread = std::thread([this, cfg, specs = std::move(specs)]() {
        try {
            const Grid& doseGrid = *m_state.computeGrid;
            auto objectives = ObjectiveBuilder::buildFromSpecs(
                specs, *m_state.patientData, m_state.patientData->getGrid(), doseGrid);
            m_pending = std::make_shared<ArcPipelineResult>(ArcPipeline::run(
                *m_state.stf, *m_state.plan, *m_state.patientData, doseGrid, std::move(objectives), cfg,
                [this](const std::string& stage) {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    m_stage = stage;
                },
                [this](const ArcProgressSnapshot& snapshot) {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    m_progress.push_back(snapshot);
                }));
            const auto& seq = m_pending->sequence;
            char buf[256];
            snprintf(buf, sizeof(buf), "%zu control points | %.0f MU | %.0f s delivery",
                     seq.controlPoints.size(), seq.totalMU, seq.totalTimeSec);
            m_status = buf;
        } catch (const std::exception& e) {
            Logger::error("DAO failed: " + std::string(e.what()));
            m_status = "Error: " + std::string(e.what());
            m_pending.reset();
        }
        m_done = true;
    });
}

void DaoOptimizationView::renderStatus() {
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
            if (m_pending) {
                const int number = m_state.nextLeafSeqNumber();
                const int doseId = m_state.doseStore.addDose(
                    "VMAT DAO #" + std::to_string(number), m_pending->dose, m_state.computeGrid);
                m_state.incrementLeafSeqCount();
                m_state.arcCache[doseId] = {m_pending, -1};
                m_state.syncSelectedDose();
                m_state.optimizationJustFinished = true;
                m_pending.reset();
            }
        }
        return;
    }
    if (!m_status.empty()) {
        ImGui::TextColored(m_state.arcResult ? getThemeColors().passText : getThemeColors().failText,
                           "%s", m_status.c_str());
        if (m_state.arcResult) ImGui::TextDisabled("Apertures and control points: Leaf Sequencing panel.");
    }
}

void DaoOptimizationView::renderProgress() {
    const std::vector<ArcProgressSnapshot>* historySource = nullptr;
    if (m_running.load()) {
        std::vector<ArcProgressSnapshot> updates;
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            updates.swap(m_progress);
        }
        for (auto& sample : updates) m_displayHistory.push_back(std::move(sample));
        historySource = &m_displayHistory;
    } else if (m_state.arcResult) {
        historySource = &m_state.arcResult->progressHistory;
    }
    if (!historySource || historySource->empty()) return;
    const auto& history = *historySource;
    if (!ImGui::CollapsingHeader("DAO Progress", ImGuiTreeNodeFlags_DefaultOpen)) return;

    const auto& latest = history.back();
    m_progressStage = std::min(m_progressStage, latest.stage);
    const std::string selection = m_progressStage == 0 ? "Latest resolution" :
        "Resolution " + std::to_string(m_progressStage);
    if (ImGui::BeginCombo("History", selection.c_str())) {
        if (ImGui::Selectable("Latest resolution", m_progressStage == 0)) m_progressStage = 0;
        for (int stage = 1; stage <= latest.stage; ++stage) {
            const auto label = "Resolution " + std::to_string(stage);
            if (ImGui::Selectable(label.c_str(), m_progressStage == stage)) m_progressStage = stage;
        }
        ImGui::EndCombo();
    }
    const int stage = m_progressStage == 0 ? latest.stage : std::min(m_progressStage, latest.stage);
    ImGui::Checkbox("Log cost", &m_logCost);
    if (m_logCost) ImGui::TextDisabled("log10(cost), floored at 1e-12 for zero costs.");
    std::vector<float> costs, mu;
    const ArcProgressSnapshot* selected = nullptr;
    const ArcProgressSnapshot* lastStep = nullptr;
    bool invalid = false;
    for (const auto& snapshot : history) {
        if (snapshot.stage != stage) continue;
        selected = &snapshot;
        if (!snapshot.stageFinished && snapshot.progress.acceptedSteps > 0) lastStep = &snapshot;
        const double cost = snapshot.breakdown.total;
        if (!std::isfinite(cost) || !std::isfinite(snapshot.progress.totalMu) || cost < 0.0 ||
            snapshot.progress.totalMu < 0.0 || (!m_logCost && cost > FLT_MAX) ||
            snapshot.progress.totalMu > FLT_MAX) {
            invalid = true;
            continue;
        }
        costs.push_back(static_cast<float>(m_logCost ? std::log10(std::max(cost, 1e-12)) : cost));
        mu.push_back(static_cast<float>(snapshot.progress.totalMu));
    }
    if (!selected) return;
    const auto& sample = *selected;
    ImGui::Text("Resolution %d/%d | %zu control points | iteration %d/%d",
        sample.stage, sample.stages, sample.controlPoints, sample.progress.iteration, sample.maxIterations);
    const float fraction = sample.stageFinished ? 1.0f :
        static_cast<float>(sample.progress.iteration) / std::max(sample.maxIterations, 1);
    ImGui::ProgressBar(std::clamp(fraction, 0.0f, 1.0f), ImVec2(-1, 0),
                      sample.stageFinished ? "Resolution finished" : "Resolution iterations");
    ImGui::Text("Cost %.4e | MU %.1f | max dose %.2f Gy",
                sample.breakdown.total, sample.progress.totalMu, sample.breakdown.maxDose);
    ImGui::Text("Accepted %d | elapsed %.1f s", sample.progress.acceptedSteps, sample.progress.elapsedSec);
    if (lastStep && ImGui::TreeNode("Last sampled step")) {
        const auto& info = lastStep->progress;
        ImGui::Text("Iteration %d | accepted %d", info.iteration, info.acceptedSteps);
        ImGui::Text("GradInf %.4e", info.preStepGradInf);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Pre-step gradient infinity norm in scaled leaf/MU variables; not projected.");
        ImGui::Text("Step %.2e | LS iterations %d", info.stepSize, info.lsIters);
        ImGui::TreePop();
    }
    ImGui::TextDisabled("Sampled updates: initial, first accepted step, every 10 steps, final.");
    if (invalid) {
        ImGui::TextColored(getThemeColors().failText, "Invalid cost/MU samples; charts unavailable.");
    } else if (!costs.empty()) {
        ImGui::PlotLines("Objective##dao", costs.data(), static_cast<int>(costs.size()), 0,
                         nullptr, FLT_MAX, FLT_MAX, ImVec2(0, 130));
        ImGui::PlotLines("Total MU##dao", mu.data(), static_cast<int>(mu.size()), 0,
                         nullptr, FLT_MAX, FLT_MAX, ImVec2(0, 100));
        ImGui::TextDisabled("X axis: sampled updates within this resolution (not uniform iterations).");
    }

    ImGui::SeparatorText("Weighted penalty contributions");
    ImGui::TextDisabled("Cost share, not clinical goal attainment or gradient magnitude.");
    if (!std::isfinite(sample.breakdown.total) || sample.breakdown.total < 0.0) {
        ImGui::TextColored(getThemeColors().failText, "Invalid total cost; bars unavailable.");
        return;
    }
    auto terms = sample.breakdown.contributions;
    for (const auto& term : terms) {
        if (!std::isfinite(term.second) || term.second < 0.0) {
            ImGui::TextColored(getThemeColors().failText, "Invalid penalty value; bars unavailable.");
            return;
        }
    }
    std::stable_sort(terms.begin(), terms.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    for (size_t i = 0; i < terms.size(); ++i) {
        ImGui::PushID(static_cast<int>(i));
        ImGui::TextWrapped("%s", terms[i].first.c_str());
        const double share = sample.breakdown.total > 0.0 ? terms[i].second / sample.breakdown.total : 0.0;
        char label[80];
        snprintf(label, sizeof(label), "%.2f%% | cost %.4e", 100.0 * share, terms[i].second);
        if (!std::isfinite(share) || !std::isfinite(terms[i].second) || terms[i].second < 0.0)
            ImGui::TextColored(getThemeColors().failText, "Invalid penalty value");
        else
            ImGui::ProgressBar(static_cast<float>(std::clamp(share, 0.0, 1.0)), ImVec2(-1, 0), label);
        ImGui::PopID();
    }
}

} // namespace optirad
