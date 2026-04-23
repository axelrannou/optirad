#pragma once

#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/Machine.hpp"
#include "steering/StfProperties.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "dose/DoseStore.hpp"
#include "core/workflow/PlanAnalysis.hpp"
#include "geometry/Grid.hpp"
#include "core/Aperture.hpp"
#include <memory>
#include <vector>
#include <string>

namespace optirad {

/// Shared pipeline state used by both CLI and GUI.
/// Holds the data flowing through the treatment planning workflow:
///   load DICOM → create Plan → generate STF → dose calc → optimize.
struct WorkflowState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;
    std::shared_ptr<Stf> stf;

    /// One PhaseSpaceBeamSource per beam (one per gantry angle in the plan)
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;

    // ── Dose calculation state ──
    std::shared_ptr<DoseInfluenceMatrix> dij;
    std::shared_ptr<Grid> displayGrid;  // Grid synced with doseResult for display/analysis
    std::shared_ptr<Grid> computeGrid;

    // ── Optimization state ──
    std::vector<double> optimizedWeights;
    std::shared_ptr<DoseMatrix> doseResult;

    // ── Leaf sequencing state ──
    std::vector<LeafSequenceResult> leafSequences;
    std::vector<double> deliverableWeights;
    std::shared_ptr<DoseMatrix> deliverableDoseResult;
    std::vector<StructureDoseStats> deliverableStats;

    // ── Multi-dose management ──
    DoseStore doseStore;

    // ── Workflow queries ──
    bool dicomLoaded() const { return patientData != nullptr; }
    bool planCreated() const { return plan != nullptr; }
    bool stfGenerated() const { return stf != nullptr && !stf->isEmpty(); }
    bool phaseSpaceLoaded() const {
        return !phaseSpaceSources.empty() && phaseSpaceSources[0]->isBuilt();
    }
    bool dijComputed() const { return dij != nullptr && dij->getNumNonZeros() > 0; }
    bool optimizationDone() const { return !optimizedWeights.empty() && doseResult != nullptr; }
    bool leafSequenceDone() const { return !leafSequences.empty(); }
    bool doseAvailable() const { return doseStore.count() > 0 && doseStore.getSelected() != nullptr; }

    bool isPhaseSpaceMachine() const {
        return plan && plan->getMachine().isPhaseSpace();
    }

    /// Sync doseResult/displayGrid from DoseStore's current selection.
    void syncSelectedDose() {
        auto* sel = doseStore.getSelected();
        if (sel) {
            doseResult = sel->dose;
            displayGrid = sel->grid;
        } else {
            doseResult.reset();
        }
    }

    // ── Reset cascades ──
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
        syncSelectedDose();
    }
    void resetLeafSequence() {
        leafSequences.clear();
        deliverableWeights.clear();
        deliverableDoseResult.reset();
        deliverableStats.clear();
    }
    void resetOptimization() {
        optimizedWeights.clear();
        resetLeafSequence();
        syncSelectedDose();
    }
    void resetAllDoses() {
        doseStore.clear();
        doseResult.reset();
        displayGrid.reset();
        computeGrid.reset();
    }
};

} // namespace optirad
