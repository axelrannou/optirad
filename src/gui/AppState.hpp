#pragma once

#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/Machine.hpp"
#include "steering/StfProperties.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include "DoseManager.hpp"
#include <memory>
#include <vector>
#include <string>
#include <atomic>

namespace optirad {

/// Shared application state for the GUI pipeline.
/// Mirrors the CLI AppState: load DICOM → create Plan → generate STF → dose calc → optimize.
struct GuiAppState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;
    std::shared_ptr<Stf> stf;

    /// One PhaseSpaceBeamSource per beam (one per gantry angle in the plan)
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;

    // ── Dose calculation state ──
    std::shared_ptr<DoseInfluenceMatrix> dij;
    std::shared_ptr<Grid> doseGrid;        // Display grid (synced with doseResult for rendering)
    std::shared_ptr<Grid> computeGrid;     // Computation grid for Dij / optimization

    // ── Optimization state ──
    std::vector<double> optimizedWeights;
    std::shared_ptr<DoseMatrix> doseResult;

    // ── Multi-dose management ──
    DoseManager doseManager;

    // ── Async task state ──
    std::atomic<bool> cancelFlag{false};
    std::string taskStatus;   // e.g. "Calculating Dij..."
    float taskProgress = 0.f; // 0..1
    bool taskRunning = false;

    /// Set to true when optimization finishes; Application reads & resets it to switch DVH tab.
    std::atomic<bool> optimizationJustFinished{false};

    // Workflow state
    bool dicomLoaded() const { return patientData != nullptr; }
    bool planCreated() const { return plan != nullptr; }
    bool stfGenerated() const { return stf != nullptr && !stf->isEmpty(); }
    bool phaseSpaceLoaded() const {
        return !phaseSpaceSources.empty() && phaseSpaceSources[0]->isBuilt();
    }
    bool dijComputed() const { return dij != nullptr && dij->getNumNonZeros() > 0; }
    bool optimizationDone() const { return !optimizedWeights.empty() && doseResult != nullptr; }
    bool doseAvailable() const { return doseManager.count() > 0 && doseManager.getSelected() != nullptr; }

    /// Check if the current plan uses a phase-space machine
    bool isPhaseSpaceMachine() const {
        return plan && plan->getMachine().isPhaseSpace();
    }

    /// Sync doseResult/doseGrid from DoseManager's current selection (for backward compat).
    void syncSelectedDose() {
        auto* sel = doseManager.getSelected();
        if (sel) {
            doseResult = sel->dose;
            doseGrid = sel->grid;
        } else {
            doseResult.reset();
            // Don't clear doseGrid here — it may still be used by the dose engine
        }
    }

    // Reset downstream state when upstream changes.
    // After clearing pipeline state, re-sync display dose from DoseManager so
    // SliceViews / stats / DVH keep showing the currently selected dose.
    void resetPlan() {
        plan.reset(); stfProps.reset(); stf.reset(); phaseSpaceSources.clear();
        resetDij();
    }
    void resetStf() { stfProps.reset(); stf.reset(); resetDij(); }
    void resetPhaseSpace() { phaseSpaceSources.clear(); }
    void resetDij() {
        dij.reset();
        computeGrid.reset();
        optimizedWeights.clear();
        // Restore display dose from DoseManager.
        syncSelectedDose();
    }
    void resetOptimization() {
        optimizedWeights.clear();
        syncSelectedDose();
    }
    void resetAllDoses() { doseManager.clear(); doseResult.reset(); doseGrid.reset(); computeGrid.reset(); }
};

} // namespace optirad
