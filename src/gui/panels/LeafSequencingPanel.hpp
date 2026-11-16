#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
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
    /// VMAT plans: apertures and control points of the selected DAO result (DAO runs in the Optimization panel).
    void renderArc();

    GuiAppState& m_state;

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
