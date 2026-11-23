#pragma once

#include "DeliveryPlayback.hpp"

#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/Machine.hpp"
#include "steering/StfProperties.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/IDoseOperator.hpp"
#include "dose/DoseMatrix.hpp"
#include "dose/DoseStore.hpp"
#include "core/workflow/PlanAnalysis.hpp"
#include "geometry/Grid.hpp"
#include "core/Aperture.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "core/workflow/DoseCalculationPipeline.hpp"
#include <array>
#include "optimization/ObjectiveProtocol.hpp"
#include "io/DicomContext.hpp"
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

/// Cached arc (VMAT) results, keyed by the arc dose entry ID.
struct ArcCacheEntry {
    std::shared_ptr<ArcPipelineResult> result;
    int linkedOptDoseId = -1; // the fluence optimization dose it was derived from
};

/// Shared application state for the GUI pipeline.
/// Mirrors the CLI AppState: load DICOM → create Plan → generate STF → dose calc → optimize.
struct GuiAppState {
    DeliveryPlayback delivery;
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;
    std::shared_ptr<Stf> stf;

    /// One PhaseSpaceBeamSource per beam (one per gantry angle in the plan)
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;

    // ── Dose calculation state ──
    std::shared_ptr<DoseInfluenceMatrix> dij;       // stored-Dij engine only
    std::shared_ptr<IDoseOperator> doseOp;          // dij or matrix-free operator
    DoseComputeBackend computeBackend = DoseComputeBackend::Auto;
    size_t gpuMemoryBudgetMiB = 0;
    std::shared_ptr<Grid> displayGrid;      // Grid synced with doseResult for rendering
    std::shared_ptr<Grid> computeGrid;      // Computation grid for Dij / optimization

    /// Dose settings of a DAO plan, which needs no Dij or fluence model: only the dose grid and the engine settings.
    struct DaoDoseSettings {
        std::array<double, 3> resolution = DoseCalcPipelineOptions{}.resolution;
        int numThreads = 0;
        bool excludeExternal = DoseCalcPipelineOptions{}.excludeExternal;
        double fluencePixelMm = DoseCalcPipelineOptions{}.fluencePixelMm;
        DoseComputeBackend backend = DoseComputeBackend::Auto;
        size_t gpuMemoryBudgetMiB = 0;
    } daoDose;

    // ── Optimization state ──
    std::vector<double> optimizedWeights;
    std::shared_ptr<DoseMatrix> doseResult;

    // ── Leaf sequencing state ──
    std::vector<LeafSequenceResult> leafSequences;
    std::vector<double> deliverableWeights;
    std::shared_ptr<DoseMatrix> deliverableDoseResult;
    std::vector<StructureDoseStats> deliverableStats;

    // ── Arc (VMAT) state ──
    std::shared_ptr<ArcPipelineResult> arcResult;
    std::vector<ObjectiveSpec> optimizationSpecs;   // objectives of the last fluence optimization (reused by the arc stage)
    OptimizationConfig optimizationConfig;          // settings of the last fluence optimization

    // ── Multi-dose management ──
    DoseStore doseStore;

    // ── DICOM context (captured during import for use in export) ──
    DicomContext dicomContext;

    // ── Output directory shared by Export Dose / RT Plan / Dose Stats CSV ──
    std::string exportDir;

    // ── Pipeline result caches (keyed by dose entry ID) ──
    std::unordered_map<int, std::vector<double>> optWeightsCache;
    std::unordered_map<int, LeafSeqCacheEntry> seqCache;
    std::unordered_map<int, ArcCacheEntry> arcCache;
    std::unordered_map<int, std::vector<StructureDoseStats>> statsCache;
    std::vector<int> statsDoseIds;
    const PatientData* statsPatient = nullptr;
    double statsPrescribedDose = 0.0;
    int activeOptDoseId = -1; // which optimization dose produced current weights

    // ── Plan of each dose ──
    // Every dose is delivered by a plan (null: none, e.g. an RTDOSE imported without its RTPLAN). Selecting a
    // dose makes its plan the current one, so the views, the playback and the leaf panel show that plan.
    struct PlanContext {
        std::shared_ptr<Plan> plan;
        std::shared_ptr<StfProperties> stfProps;
        std::shared_ptr<Stf> stf;
        std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;
        std::vector<ObjectiveSpec> optimizationSpecs;
        std::shared_ptr<DoseInfluenceMatrix> dij;   // kept for the working plan only
        std::shared_ptr<IDoseOperator> doseOp;
        std::shared_ptr<Grid> computeGrid;
    };
    std::unordered_map<int, std::shared_ptr<Plan>> dosePlans;     // dose entry ID -> plan that delivers it
    std::unordered_map<const Plan*, PlanContext> parkedPlans;      // plans that are not the current one
    std::shared_ptr<Plan> workingPlan;     // last plan created in the Planning panel (or imported)
    std::shared_ptr<Plan> importedPlan;    // plan built from the imported RTPLAN
    std::shared_ptr<ArcPipelineResult> importedDelivery;  // its leaves and MU (arc RTPLAN), without dose
    std::vector<LeafSequenceResult> importedSequences;     // its segments and MU (static-beam RTPLAN)

    // ── Workflow naming counters ──
    int optimizationCount = 0;
    int leafSeqCount = 0;

    int nextOptimizationNumber() const { return optimizationCount + 1; }
    void incrementOptimizationCount() { ++optimizationCount; }
    int nextLeafSeqNumber() const { return leafSeqCount + 1; }
    void incrementLeafSeqCount() { ++leafSeqCount; }

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
    bool dijComputed() const { return doseOp != nullptr && (!dij || dij->getNumNonZeros() > 0); }
    /// DAO on a VMAT plan skips the Dij / fluence optimization stages.
    bool directDao() const {
        return plan && plan->isVmat();
    }
    bool optimizationDone() const { return !optimizedWeights.empty() && doseResult != nullptr; }
    bool leafSequenceDone() const { return !leafSequences.empty(); }
    /// Static-beam leaf sequences to display: the selected dose's sequencing result, or the imported RTPLAN's
    /// segments while its plan is current (also without an imported dose).
    const std::vector<LeafSequenceResult>* deliveredSequences() const {
        if (const auto* sel = doseStore.getSelected()) {
            const auto it = seqCache.find(sel->id);
            if (it != seqCache.end()) return &it->second.sequences;
        }
        if (plan && plan == importedPlan && !importedSequences.empty()) return &importedSequences;
        return nullptr;
    }
    bool doseAvailable() const { return doseStore.count() > 0 && doseStore.getSelected() != nullptr; }

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

    /// Analyse every loaded dose on the same CT samples in their common grid coverage.
    const std::vector<StructureDoseStats>& getOrComputeStats(
            int entryIdx, const PatientData& patient, double prescribedDose = 60.0) {
        static const std::vector<StructureDoseStats> empty;
        auto* entry = doseStore.getEntry(entryIdx);
        if (!entry || !entry->dose || !entry->grid) return empty;
        std::vector<int> doseIds;
        std::vector<Grid> comparisonGrids;
        for (const auto& candidate : doseStore.getEntries()) {
            if (!candidate.dose || !candidate.grid) continue;
            doseIds.push_back(candidate.id);
            comparisonGrids.push_back(*candidate.grid);
        }
        if (doseIds != statsDoseIds || statsPatient != &patient ||
            statsPrescribedDose != prescribedDose) {
            statsCache.clear();
            statsDoseIds = std::move(doseIds);
            statsPatient = &patient;
            statsPrescribedDose = prescribedDose;
        }
        auto it = statsCache.find(entry->id);
        if (it != statsCache.end()) return it->second;
        auto stats = PlanAnalysis::computeStats(
            *entry->dose, patient, *entry->grid, prescribedDose, comparisonGrids);
        auto [ins, _] = statsCache.emplace(entry->id, std::move(stats));
        return ins->second;
    }

    void invalidateStatsCache(int entryId) { statsCache.erase(entryId); }

    /// Sync doseResult/displayGrid from DoseStore's current selection,
    /// and restore cached optimization / leaf sequencing results.
    void syncSelectedDose() {
        auto* sel = doseStore.getSelected();
        if (sel) {
            // A new dose is delivered by the current plan (doses are added right after they are computed).
            auto linked = dosePlans.try_emplace(sel->id, plan).first->second;
            switchToPlan(linked ? linked : (workingPlan ? workingPlan : plan));
        } else if (workingPlan) {
            switchToPlan(workingPlan);
        }
        if (sel) {
            doseResult = sel->dose;
            displayGrid = sel->grid;
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
            arcResult.reset();
            for (const auto& [arcDoseId, arcEntry] : arcCache) {
                if (arcEntry.linkedOptDoseId == selId) { arcResult = arcEntry.result; break; }
            }
            return;
        }

        // Selected dose is an arc dose: restore its result and the linked optimization weights
        auto arcIt = arcCache.find(selId);
        if (arcIt != arcCache.end()) {
            arcResult = arcIt->second.result;
            auto linkedIt = optWeightsCache.find(arcIt->second.linkedOptDoseId);
            if (linkedIt != optWeightsCache.end()) {
                optimizedWeights = linkedIt->second;
                activeOptDoseId = arcIt->second.linkedOptDoseId;
            } else {
                optimizedWeights.clear();
                activeOptDoseId = -1;
            }
            leafSequences.clear();
            deliverableWeights.clear();
            deliverableStats.clear();
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

        // No cached results for this dose — clear pipeline state (an imported RTPLAN keeps its delivery)
        const bool importedCurrent = plan && plan == importedPlan;
        arcResult = importedCurrent ? importedDelivery : nullptr;
        optimizedWeights.clear();
        if (importedCurrent) leafSequences = importedSequences;
        else leafSequences.clear();
        deliverableWeights.clear();
        deliverableStats.clear();
        activeOptDoseId = -1;
    }

    bool planHasDose(const Plan* p) const {
        for (const auto& [id, linked] : dosePlans)
            if (linked.get() == p) return true;
        return false;
    }

    /// Make `target` the current plan: the current plan is parked (if a dose or the working plan uses it) and
    /// the target's parked state (beams, objectives, Dij of the working plan) restored.
    void switchToPlan(const std::shared_ptr<Plan>& target) {
        if (target == plan) return;
        if (plan && (plan == workingPlan || plan == importedPlan || planHasDose(plan.get()))) {
            PlanContext ctx{plan, stfProps, stf, phaseSpaceSources, optimizationSpecs, nullptr, nullptr, nullptr};
            if (plan == workingPlan) { ctx.dij = dij; ctx.doseOp = doseOp; ctx.computeGrid = computeGrid; }
            parkedPlans[plan.get()] = std::move(ctx);
        }
        PlanContext next;
        next.plan = target;
        if (auto it = parkedPlans.find(target.get()); target && it != parkedPlans.end()) {
            next = std::move(it->second);
            parkedPlans.erase(it);
        }
        plan = next.plan;
        stfProps = next.stfProps;
        stf = next.stf;
        phaseSpaceSources = std::move(next.phaseSpaceSources);
        optimizationSpecs = std::move(next.optimizationSpecs);
        dij = next.dij;
        doseOp = next.doseOp;
        computeGrid = next.computeGrid;
        optimizedWeights.clear();
    }

    /// Clear the current plan before a new one is created: no dose is selected, so the new plan starts without
    /// dose; the plans of existing doses stay with them.
    void startNewPlan() {
        workingPlan.reset();
        doseStore.clearSelection();
        switchToPlan(nullptr);
        resetLeafSequence();
        syncSelectedDose();
    }

    // Reset downstream state when upstream changes.
    // After clearing pipeline state, re-sync display dose from DoseStore so
    // SliceViews / stats / DVH keep showing the currently selected dose.
    void resetPlan() {
        plan.reset(); stfProps.reset(); stf.reset(); phaseSpaceSources.clear(); optimizationSpecs.clear();
        resetDij();
    }
    void resetStf() { stfProps.reset(); stf.reset(); resetDij(); }
    void resetPhaseSpace() { phaseSpaceSources.clear(); }
    void resetDij() {
        dij.reset();
        doseOp.reset();
        computeGrid.reset();
        optimizedWeights.clear();
        // Restore display dose from DoseStore.
        syncSelectedDose();
    }
    void resetLeafSequence() {
        arcResult.reset();
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
        doseStore.clear(); doseResult.reset(); displayGrid.reset(); computeGrid.reset();
        optWeightsCache.clear(); seqCache.clear(); arcCache.clear(); statsCache.clear(); activeOptDoseId = -1;
        dosePlans.clear(); parkedPlans.clear(); workingPlan.reset(); importedPlan.reset(); importedDelivery.reset();
        importedSequences.clear();
        statsDoseIds.clear(); statsPatient = nullptr; statsPrescribedDose = 0.0;
    }
};

} // namespace optirad
