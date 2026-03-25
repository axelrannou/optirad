#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/workflow/OptimizationPipeline.hpp"
#include "optimization/ObjectiveProtocol.hpp"
#include <string>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>

namespace optirad {

/// Optimization panel: configure objectives per structure, run L-BFGS-B optimizer.
/// Activates once Dij is computed.
class OptimizationPanel : public IPanel {
public:
    OptimizationPanel(GuiAppState& state);
    ~OptimizationPanel();

    std::string getName() const override { return "Optimization"; }
    void update() override {}
    void render() override;

private:
    GuiAppState& m_state;

    // Per-structure objective configuration
    struct ObjectiveConfig {
        std::string structureName;
        int typeIdx = 0; // 0=SquaredDeviation, 1=SquaredOverdose, 2=SquaredUnderdose, 3=MinDVH, 4=MaxDVH
        float doseValue = 60.0f;     // Gy
        float weight = 1.0f;
        float volumePct = 95.0f;     // Volume fraction (%) for DVH objectives
    };
    std::vector<ObjectiveConfig> m_objectives;
    bool m_objectivesInitialized = false;

    // Optimizer settings
    int m_maxIterations = 400;
    float m_tolerance = 1e-5f;

    // NTO (Normal Tissue Objective) / hotspot control
    bool m_ntoEnabled = true;
    float m_ntoThresholdPct = 104.0f;  // % of prescription dose
    int m_ntoPenalty = 2000.0f;

    // Async state
    std::atomic<bool> m_isOptimizing{false};
    std::atomic<bool> m_optimizationDone{false};
    std::thread m_optThread;
    std::string m_optStatusMessage;
    std::atomic<int> m_currentIteration{0};
    std::atomic<double> m_currentObjective{0.0};
    std::atomic<double> m_currentProjGrad{0.0};
    std::atomic<double> m_currentImprovement{0.0};
    OptimizationPipelineResult m_pipelineResult;

    // Iteration log for convergence curve (written by optimizer thread, read by render)
    std::mutex m_iterMutex;
    std::vector<IterationInfo> m_iterationLog;

    void renderConvergenceCurve();
};

} // namespace optirad
