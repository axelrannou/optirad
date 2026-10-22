#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "core/workflow/LeafSequencingPipeline.hpp"
#include <string>
#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>

namespace optirad {

/// Leaf sequencing panel: configure intensity levels, run step-and-shoot
/// decomposition, display per-beam results and deliverable dose statistics.
/// Activates once optimization is complete.
class LeafSequencingPanel : public IPanel {
public:
    LeafSequencingPanel(GuiAppState& state);
    ~LeafSequencingPanel();

    std::string getName() const override { return "Leaf Sequencing"; }
    void update() override {}
    void render() override;

private:
    /// VMAT plans: arc sequencing + aperture weight optimization instead of step-and-shoot sequencing.
    void renderArc();

    GuiAppState& m_state;

    // VMAT arc stage (initial values come from the pipeline defaults, shared with the CLI)
    int m_arcStrategy = static_cast<int>(ArcOptimizationConfig{}.strategy);   // SweepRefine, ProgressiveDao
    int m_arcDaoIterations = ArcOptimizationConfig{}.daoIterations;
    float m_arcDaoCoarsest = static_cast<float>(ArcOptimizationConfig{}.daoCoarsestSpacingDeg);
    float m_arcDaoHotspotThresholdPct = static_cast<float>(100.0 * ArcOptimizationConfig{}.daoHotspotThresholdPct);
    int m_arcDaoHotspotPenalty = static_cast<int>(ArcOptimizationConfig{}.daoHotspotPenalty);
    float m_arcDaoMuPenalty = static_cast<float>(ArcOptimizationConfig{}.daoMuPenalty);
    int m_arcMethod = static_cast<int>(ArcOptimizationConfig{}.sequencer.method);   // Auto, Sweep, LevelSet
    int m_arcLevels = ArcOptimizationConfig{}.sequencer.numLevels;
    float m_arcRiseRatio = static_cast<float>(ArcOptimizationConfig{}.sequencer.maxRiseRatio);
    float m_arcSweepFraction = static_cast<float>(ArcOptimizationConfig{}.sequencer.maxSweepFraction);
    bool m_arcIntegrateMotion = ArcOptimizationConfig{}.integrateLeafMotion;
    bool m_arcOptimizeWeights = ArcOptimizationConfig{}.optimizeWeights;
    bool m_arcOptimizeLeaves = ArcOptimizationConfig{}.optimizeLeaves;
    int m_arcLeafIterations = ArcOptimizationConfig{}.leafIterations;
    int m_arcIterations = ArcOptimizationConfig{}.weightIterations;
    float m_arcMuSmoothing = static_cast<float>(ArcOptimizationConfig{}.muSmoothingWeight);
    float m_arcHotspotThresholdPct = static_cast<float>(100.0 * ArcOptimizationConfig{}.hotspotThresholdPct);
    int m_arcHotspotPenalty = static_cast<int>(ArcOptimizationConfig{}.hotspotPenalty);
    float m_arcMuPenalty = static_cast<float>(ArcOptimizationConfig{}.muPenalty);
    float m_arcMaxWeightFactor = static_cast<float>(ArcOptimizationConfig{}.maxWeightFactor);
    std::atomic<bool> m_arcRunning{false};
    std::atomic<bool> m_arcDone{false};
    std::thread m_arcThread;
    std::mutex m_arcMutex;           // guards m_arcStage
    std::string m_arcStage;
    std::atomic<int> m_arcIter{0};
    std::shared_ptr<ArcPipelineResult> m_arcPending;
    std::string m_arcStatus;

    // Sequencer settings
    int m_numLevels = LeafSequencerOptions{}.numLevels;
    float m_minSegmentMU = static_cast<float>(LeafSequencerOptions{}.minSegmentMU);

    // Async state
    std::atomic<bool> m_isRunning{false};
    std::atomic<bool> m_isDone{false};
    std::thread m_thread;
    std::string m_statusMessage;
    std::atomic<int> m_currentBeam{0};
    std::atomic<int> m_totalBeams{0};
    LeafSequencingPipelineResult m_pipelineResult;

    // Display cache: status message per deliverable dose ID
    struct SeqDisplayCache {
        std::string statusMessage;
    };
    std::unordered_map<int, SeqDisplayCache> m_displayCache;
    int m_doseVersion = -1;
};

} // namespace optirad
