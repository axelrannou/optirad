#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/workflow/LeafSequencingPipeline.hpp"
#include <string>
#include <atomic>
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
    GuiAppState& m_state;

    // Sequencer settings
    int m_numLevels = 15;
    float m_minSegmentMU = 0.0f;

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
