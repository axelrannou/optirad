#pragma once

#include "optimization/ObjectiveBuilder.hpp"
#include "optimization/ObjectiveProtocol.hpp"
#include "optimization/IOptimizer.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "dose/PlanAnalysis.hpp"
#include "core/PatientData.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <vector>
#include <string>
#include <atomic>
#include <functional>

namespace optirad {

/// Configuration for the optimization pipeline.
struct OptimizationConfig {
    int maxIterations = 400;
    double tolerance = 1e-5;
    double targetDose = 66.0;

    // NTO / hotspot control
    bool ntoEnabled = true;
    double ntoThresholdPct = 1.04;  // fraction of Rx (104%)
    double ntoPenalty = 2000.0;
};

/// Result of the optimization pipeline.
struct OptimizationPipelineResult {
    std::vector<double> weights;
    std::shared_ptr<DoseMatrix> doseResult;
    std::vector<StructureDoseStats> stats;
    int iterations = 0;
    double finalObjective = 0.0;
    bool converged = false;
};

/// Extracts the shared optimization workflow from CLI optimize() and GUI OptimizationPanel.
/// Handles: objective building → optimizer → forward dose → plan analysis.
class OptimizationPipeline {
public:
    /// Run optimization using a protocol (auto-detects PTV targets).
    static OptimizationPipelineResult run(
        const DoseInfluenceMatrix& dij,
        const OptimizationConfig& config,
        const ObjectiveProtocol& protocol,
        const PatientData& patientData,
        const Grid& doseGrid);

    /// Run optimization using pre-built objectives.
    static OptimizationPipelineResult runWithObjectives(
        const DoseInfluenceMatrix& dij,
        const OptimizationConfig& config,
        BuiltObjectives objectives,
        const PatientData& patientData,
        const Grid& doseGrid);
};

} // namespace optirad
