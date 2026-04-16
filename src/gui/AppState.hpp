#pragma once

#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/Machine.hpp"
#include "steering/StfProperties.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "dose/DoseManager.hpp"
#include "dose/PlanAnalysis.hpp"
#include "geometry/Grid.hpp"
#include "core/Aperture.hpp"
#include <memory>
#include <vector>
#include <string>
#include <atomic>
#include <unordered_map>

namespace optirad {

/// Cached leaf sequencing results, keyed by deliverable dose entry ID.
struct LeafSeqCacheEntry {
    std::vector<LeafSequenceResult> sequences;
    std::vector<double> deliverableWeights;
    std::vector<StructureDoseStats> deliverableStats;
    int linkedOptDoseId = -1; // the optimization dose these were derived from
};

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

    // ── Leaf sequencing state ──
    std::vector<LeafSequenceResult> leafSequences;
    std::vector<double> deliverableWeights;
    std::shared_ptr<DoseMatrix> deliverableDoseResult;
    std::vector<StructureDoseStats> deliverableStats;

    // ── Multi-dose management ──
    DoseManager doseManager;

    // ── Pipeline result caches (keyed by dose entry ID) ──
    std::unordered_map<int, std::vector<double>> optWeightsCache;
    std::unordered_map<int, LeafSeqCacheEntry> seqCache;
    int activeOptDoseId = -1; // which optimization dose produced current weights

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
    bool leafSequenceDone() const { return !leafSequences.empty(); }
    bool doseAvailable() const { return doseManager.count() > 0 && doseManager.getSelected() != nullptr; }

    /// Check if the current plan uses a phase-space machine
    bool isPhaseSpaceMachine() const {
        return plan && plan->getMachine().isPhaseSpace();
    }

    /// Cache current optimization weights for a dose entry ID.
    void cacheOptimization(int doseId) {
        optWeightsCache[doseId] = optimizedWeights;
        activeOptDoseId = doseId;
    }

    /// Cache current leaf sequencing results for a deliverable dose entry ID.
    void cacheLeafSequencing(int doseId, int linkedOptDoseId) {
        LeafSeqCacheEntry entry;
        entry.sequences = leafSequences;
        entry.deliverableWeights = deliverableWeights;
        entry.deliverableStats = deliverableStats;
        entry.linkedOptDoseId = linkedOptDoseId;
        seqCache[doseId] = std::move(entry);
    }

    /// Sync doseResult/doseGrid from DoseManager's current selection,
    /// and restore cached optimization / leaf sequencing results.
    void syncSelectedDose() {
        auto* sel = doseManager.getSelected();
        if (sel) {
            doseResult = sel->dose;
            doseGrid = sel->grid;
        } else {
            doseResult.reset();
        }

        int selId = sel ? sel->id : -1;

        // Try to restore optimization weights for the selected dose
        auto optIt = optWeightsCache.find(selId);
        if (optIt != optWeightsCache.end()) {
            // Selected dose IS an optimization dose
            optimizedWeights = optIt->second;
            activeOptDoseId = selId;

            // Check if leaf sequencing was done for this optimization
            // (reverse-scan: find seqCache entry whose linkedOptDoseId matches)
            bool foundSeq = false;
            for (const auto& [seqDoseId, seqEntry] : seqCache) {
                if (seqEntry.linkedOptDoseId == selId) {
                    leafSequences = seqEntry.sequences;
                    deliverableWeights = seqEntry.deliverableWeights;
                    deliverableStats = seqEntry.deliverableStats;
                    foundSeq = true;
                    break;
                }
            }
            if (!foundSeq) {
                leafSequences.clear();
                deliverableWeights.clear();
                deliverableStats.clear();
            }
            return;
        }

        // Check if selected dose is a deliverable dose
        auto seqIt = seqCache.find(selId);
        if (seqIt != seqCache.end()) {
            const auto& entry = seqIt->second;
            leafSequences = entry.sequences;
            deliverableWeights = entry.deliverableWeights;
            deliverableStats = entry.deliverableStats;

            // Also restore the linked optimization weights
            auto linkedIt = optWeightsCache.find(entry.linkedOptDoseId);
            if (linkedIt != optWeightsCache.end()) {
                optimizedWeights = linkedIt->second;
                activeOptDoseId = entry.linkedOptDoseId;
            } else {
                optimizedWeights.clear();
                activeOptDoseId = -1;
            }
            return;
        }

        // No cached results for this dose — clear pipeline state
        optimizedWeights.clear();
        leafSequences.clear();
        deliverableWeights.clear();
        deliverableStats.clear();
        activeOptDoseId = -1;
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
        doseManager.clear(); doseResult.reset(); doseGrid.reset(); computeGrid.reset();
        optWeightsCache.clear(); seqCache.clear(); activeOptDoseId = -1;
    }
};

} // namespace optirad
