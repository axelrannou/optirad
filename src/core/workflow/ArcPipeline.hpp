#pragma once

#include "core/ArcSequence.hpp"
#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/workflow/OptimizationPipeline.hpp"
#include "core/workflow/PlanAnalysis.hpp"
#include "dose/DoseCalcOptions.hpp"
#include "dose/FluenceDoseModel.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include "optimization/ObjectiveBuilder.hpp"
#include "sequencing/ArcSequencer.hpp"
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace optirad {

enum class ArcStrategy {
    /// Fluence-based: sweep sequencing of the optimized fluence, aperture weight optimization, leaf refinement.
    SweepRefine,
    /// Progressive resolution direct aperture optimization (as in Eclipse's VMAT optimizer): conformal apertures at a
    /// few widely spaced control points are optimized (leaf positions and MU against the dose objectives), then
    /// control points are inserted by interpolation and the optimization continues, until the final spacing.
    ProgressiveDao
};

struct ArcOptimizationConfig {
    ArcStrategy strategy = ArcStrategy::SweepRefine;

    /// ProgressiveDao: control point spacing of the first resolution; each further resolution halves it until the
    /// plan's control point spacing is reached (e.g. 16, 8, 4, 2 deg).
    double daoCoarsestSpacingDeg = 16.0;
    /// ProgressiveDao: iterations of the leaf/MU optimizer at every resolution but the last, which gets twice as many
    /// (`leafIterations` is the sweep strategy's budget). The hot-spot term tightens from the fluence stage's setting to
    /// the arc setting over the resolutions. The MU penalty is off for this strategy by default: open conformal
    /// apertures need little MU (see `daoMuPenalty`).
    int daoIterations = 500;
    double daoMuPenalty = 0.0;
    /// ProgressiveDao: final hot-spot term. Less strict than the sweep's: with one opening per leaf pair a stricter
    /// term trades PTV coverage for the maximum (1.01 x Rx: max 69.8 Gy but D95 63.7 Gy on the lung case).
    double daoHotspotThresholdPct = 1.03;
    double daoHotspotPenalty = 5000.0;

    /// Targets, hotspot control and iteration limits shared with the fluence optimization.
    OptimizationConfig base;
    ArcSequencerOptions sequencer;

    /// Aperture MU may grow up to this multiple of its sequenced value.
    double maxWeightFactor = 3.0;

    /// Penalty on MU differences between neighbouring control points (smooth dose rate / gantry speed).
    double muSmoothingWeight = 0.05;

    /// Hot-spot term of the arc stages (active when `base.ntoEnabled`). Stricter than the fluence stage: the
    /// delivered dose differs from the optimized fluence dose and the difference shows up as hot spots.
    double hotspotThresholdPct = 1.01;   // fraction of the prescription
    double hotspotPenalty = 50000.0;

    /// Iterations of the aperture weight (MU) optimization. `base` supplies targets, hotspot control and tolerance.
    int weightIterations = 100;

    /// Refine the MU of each control point against the objectives. Off = deliver the sequenced MU as is.
    bool optimizeWeights = true;
    bool optimizeLeaves = true;
    /// Penalty on the total MU, as the fraction of the initial objective (of each refinement stage) that the
    /// whole sequenced MU costs. 0 = none. The dose objectives decide how much MU a plan is worth.
    double muPenalty = 0.3;
    bool integrateLeafMotion = true;   // dose of a control point = average over the leaf motion in its slice
    int leafIterations = 300;
    double leafStepMm = 1.0;   // first step and base of the per-iteration move limit (4x) of the leaf refinement
    double fluencePixelMm = 1.5;
    DoseComputeBackend computeBackend = DoseComputeBackend::Auto;
    size_t gpuMemoryBudgetMiB = 0;

    /// Dose engine settings; must match the engine that produced the fluence (IDoseOperator::engineOptions()).
    DoseCalcOptions engineOptions;
};

struct ArcPipelineResult {
    ArcSequenceResult sequence;                 // final apertures, MU and timing
    std::shared_ptr<Stf> controlPointStf;       // one beam per control point
    std::shared_ptr<DoseMatrix> dose;           // dose of the delivered arc, whole phantom/patient
    std::vector<StructureDoseStats> stats;

    int iterations = 0;
    double finalObjective = 0.0;
    bool converged = false;
    int leafIterations = 0;
    double initialLeafObjective = 0.0;
    double finalLeafObjective = 0.0;
};

/// VMAT workflow after fluence optimization: arc sequencing, then aperture weight optimization.
/// Refines aperture weights and leaf positions before recomputing the final dose.
class ArcPipeline {
public:
    /// Build one beam per control point at the given angles, with the same isocenter, bixel width and
    /// target margins as the plan's fluence-optimization beams.
    static std::shared_ptr<Stf> buildControlPointStf(const Plan& plan,
                                                     const PatientData& patientData,
                                                     const std::vector<ArcControlPoint>& controlPoints);

    /// Per-bixel weights (open fraction x MU) of every control point, indexed like `cpStf`'s global bixels.
    /// With `integrateLeafMotion` the leaves keep moving through each slice (path through the adjacent
    /// control points' tips); otherwise each control point delivers its static aperture.
    static std::vector<double> bixelWeights(const Stf& cpStf, const ArcSequenceResult& sequence,
                                            bool integrateLeafMotion = true);

    static ArcPipelineResult run(const std::vector<double>& fmoWeights,
                                 const Stf& fmoStf,
                                 const Plan& plan,
                                 const PatientData& patientData,
                                 const Grid& doseGrid,
                                 BuiltObjectives objectives,
                                 const ArcOptimizationConfig& config,
                                 IterationCallback iterCallback = nullptr,
                                 std::function<void(const std::string&)> stageCallback = nullptr);
};

} // namespace optirad
