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
#include "optimization/OptimizationLogging.hpp"
#include "sequencing/ArcSequencer.hpp"
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace optirad {

/// VMAT direct aperture optimization (DAO) settings, as in Eclipse's progressive resolution VMAT optimizer:
/// conformal apertures at a few widely spaced control points are optimized (leaf positions and MU against the dose
/// objectives), then control points are inserted by interpolation and the optimization continues, until the plan's
/// control point spacing.
struct ArcOptimizationConfig {
    /// Prescription and normal tissue objective (`base.nto`) shared with the IMRT optimizer.
    OptimizationConfig base;

    /// Control point spacing of the first resolution; each further resolution halves it until the plan's control
    /// point spacing is reached (e.g. 16, 8, 4, 2 deg).
    double coarsestSpacingDeg = 16.0;
    /// Leaf/MU iterations at every resolution but the last, which gets twice as many.
    int iterationsPerResolution = 500;
    /// Penalty on the total MU, as the fraction of the objective (at the second resolution) that the whole MU costs.
    /// 0 = none: open conformal apertures need little MU and the dose objectives decide how much MU a plan is worth.
    double muPenalty = 0.0;
    /// Penalty on MU differences between neighbouring control points (smooth dose rate / gantry speed).
    double muSmoothingWeight = 0.05;
    bool integrateLeafMotion = true;   // dose of a control point = average over the leaf motion in its slice
    double leafStepMm = 1.0;   // first step and base of the per-iteration move limit (4x)
    double fluencePixelMm = 1.5;
    /// Conformal seed apertures.
    ArcSequencerOptions seed;
    /// Backend for optimization and final delivered dose: Auto tries CUDA with logged CPU fallback;
    /// CUDA requires GPU execution; CPU disables GPU tracing and deposition.
    DoseComputeBackend computeBackend = DoseComputeBackend::Auto;
    size_t gpuMemoryBudgetMiB = 0;

    /// Dose engine settings of the arc dose.
    DoseCalcOptions engineOptions;
};

struct ArcProgressSnapshot {
    int stage = 1;
    int stages = 1;
    int maxIterations = 0;
    size_t controlPoints = 0;
    bool stageFinished = false;
    DaoProgressInfo progress;
    OptimizationBreakdown breakdown;
};

using ArcProgressCallback = std::function<void(const ArcProgressSnapshot&)>;

struct ArcPipelineResult {
    std::vector<ArcProgressSnapshot> progressHistory;
    ArcSequenceResult sequence;                 // final apertures, MU and timing
    std::shared_ptr<Stf> controlPointStf;       // one beam per control point
    std::shared_ptr<DoseMatrix> dose;           // dose of the delivered arc, whole phantom/patient
    std::vector<StructureDoseStats> stats;

    int leafIterations = 0;                     // accepted DAO steps over all resolutions
    double initialLeafObjective = 0.0;
    double finalLeafObjective = 0.0;
};

/// VMAT workflow: progressive resolution DAO from a conformal seed, then the delivered dose.
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

    /// Machine MU per fraction of one unit of control point "MU" (`Plan::machineMuPerWeight`).
    static double machineMuPerWeight(const Plan& plan);

    /// Machine constraints with the dose rate converted to weight units per second (`machineMuPerWeight`),
    /// i.e. the limits that apply to the sequence's `mu` and `doseRate`.
    static MachineConstraints deliveryConstraints(const Plan& plan);

    /// `sectorStf` has one beam per arc sector (centre angle), as built for the plan; its rays (the target
    /// projection) define the conformal seed. The normal tissue objective of `config.base.nto` is added to
    /// `objectives`.
    static ArcPipelineResult run(const Stf& sectorStf,
                                 const Plan& plan,
                                 const PatientData& patientData,
                                 const Grid& doseGrid,
                                 BuiltObjectives objectives,
                                 const ArcOptimizationConfig& config,
                                 std::function<void(const std::string&)> stageCallback = nullptr,
                                 ArcProgressCallback progressCallback = nullptr);
};

} // namespace optirad
