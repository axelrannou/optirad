#pragma once

#include "../AppState.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <thread>

namespace optirad {

/// Imported RTPLAN of the Planning panel: the plan built from it (installed at import, linked to the imported
/// RTDOSE) and the delivered dose of its own leaves and MU (no optimization), stored as a new dose linked to
/// that plan, so the Leaf Sequencing panel, the delivery playback and the BEV, 3D and slice views show it.
class ImportedPlanDoseView {
public:
    explicit ImportedPlanDoseView(GuiAppState& state);
    ~ImportedPlanDoseView();

    bool running() const { return m_running.load(); }

    /// Builds the plan of the patient's imported RTPLAN (default pencil-beam machine, RTDOSE grid) with its
    /// leaves and MU and makes it the current plan; `importedDoseId` (the imported RTDOSE, or -1) is linked to
    /// it. Does nothing without an imported RTPLAN; errors are logged and leave no plan.
    static void installImported(GuiAppState& state, int importedDoseId);

    /// Summary and "Compute dose from RTPLAN leaves". `dose` holds the dose grid and engine settings.
    void render(const GuiAppState::DaoDoseSettings& dose);

private:
    void start(const GuiAppState::DaoDoseSettings& dose);

    GuiAppState& m_state;
    bool m_useRtDoseGrid = true;

    std::atomic<bool> m_running{false};
    std::atomic<bool> m_done{false};
    std::thread m_thread;
    std::mutex m_mutex;              // guards m_stage
    std::string m_stage;
    std::string m_status;
    bool m_failed = false;

    // Worker results, installed on the GUI thread
    std::shared_ptr<Plan> m_pendingPlan;
    std::shared_ptr<ArcPipelineResult> m_pending;
    std::shared_ptr<Grid> m_pendingGrid;
};

} // namespace optirad
