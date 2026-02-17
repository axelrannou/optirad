#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include <atomic>
#include <thread>

namespace optirad {

class BeamRenderer;

/// STF generation panel: generates steering file from plan.
/// Appears only when a plan exists. "Generate STF" button triggers
/// parallel beam/ray computation. Shows progress and results.
/// Also provides per-beam visibility toggles for 3D rendering.
class StfPanel : public IPanel {
public:
    StfPanel(GuiAppState& state);

    std::string getName() const override { return "STF"; }
    void update() override {}
    void render() override;

    void setBeamRenderer(BeamRenderer* renderer) { m_beamRenderer = renderer; }

private:
    GuiAppState& m_state;
    BeamRenderer* m_beamRenderer = nullptr;

    // Generation state
    std::atomic<bool> m_isGenerating{false};
    std::atomic<bool> m_generationDone{false};
    std::thread m_generationThread;
    std::string m_statusMessage;
};

} // namespace optirad
