#pragma once

#include "../AppState.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "optimization/ObjectiveProtocol.hpp"
#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace optirad {

/// VMAT direct aperture optimization in the Optimization panel: DAO settings, the worker thread and the live
/// progress (objective / MU curves and the weighted contribution of every objective).
class DaoOptimizationView {
public:
    explicit DaoOptimizationView(GuiAppState& state);
    ~DaoOptimizationView();

    bool running() const { return m_running.load(); }

    /// DAO-specific settings (iterations, resolutions, MU terms).
    void renderSettings();
    /// Sets up the dose grid from the plan's DAO dose settings and starts the optimization.
    void start(std::vector<ObjectiveSpec> specs, const OptimizationConfig& config);
    /// Running stage, completion handling (adds the dose) and the status of the last run.
    void renderStatus();
    /// Objective / MU curves and contribution bars of the running or selected DAO result.
    void renderProgress();

private:
    GuiAppState& m_state;

    // Settings (defaults shared with the CLI through ArcOptimizationConfig)
    int m_iterations = ArcOptimizationConfig{}.iterationsPerResolution;
    float m_coarsestSpacing = static_cast<float>(ArcOptimizationConfig{}.coarsestSpacingDeg);
    float m_muPenalty = static_cast<float>(ArcOptimizationConfig{}.muPenalty);
    float m_muSmoothing = static_cast<float>(ArcOptimizationConfig{}.muSmoothing);
    bool m_integrateMotion = ArcOptimizationConfig{}.integrateLeafMotion;

    std::atomic<bool> m_running{false};
    std::atomic<bool> m_done{false};
    std::thread m_thread;
    std::mutex m_mutex;              // guards stage and live progress snapshots
    std::vector<ArcProgressSnapshot> m_progress;
    std::vector<ArcProgressSnapshot> m_displayHistory;
    std::string m_stage;
    std::shared_ptr<ArcPipelineResult> m_pending;
    std::string m_status;

    int m_progressStage = 0;         // 0 follows the latest resolution
    bool m_logCost = true;
};

} // namespace optirad
