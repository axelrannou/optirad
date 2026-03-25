#include "OptimizationPipeline.hpp"
#include "optimization/OptimizerFactory.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/Constraint.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "utils/Logger.hpp"
#include <chrono>

namespace optirad {

static OptimizationPipelineResult runImpl(
    const DoseInfluenceMatrix& dij,
    const OptimizationConfig& config,
    BuiltObjectives objectives,
    const PatientData& patientData,
    const Grid& doseGrid) {

    OptimizationPipelineResult result;

    if (objectives.ptrs.empty()) {
        Logger::error("OptimizationPipeline: no objectives to optimize");
        return result;
    }

    Logger::info("OptimizationPipeline: " + std::to_string(objectives.ptrs.size()) +
                 " objectives, bixels=" + std::to_string(dij.getNumBixels()) +
                 " voxels=" + std::to_string(dij.getNumVoxels()));

    auto optimizer = OptimizerFactory::create("LBFGS");
    optimizer->setMaxIterations(config.maxIterations);
    optimizer->setTolerance(config.tolerance);

    if (auto* lbfgs = dynamic_cast<LBFGSOptimizer*>(optimizer.get())) {
        double maxBixelDose = dij.getMaxValue();
        double estimatedMaxFluence = (maxBixelDose > 1e-10)
            ? config.targetDose / maxBixelDose * 2.0
            : 1000.0;
        lbfgs->setMaxFluence(estimatedMaxFluence);
        lbfgs->setPrescriptionDose(config.targetDose);

        if (config.ntoEnabled) {
            lbfgs->setHotspotThreshold(config.ntoThresholdPct);
            lbfgs->setHotspotPenalty(config.ntoPenalty);
        }
    }

    auto start = std::chrono::steady_clock::now();
    std::vector<Constraint> constraints;
    auto optResult = optimizer->optimize(dij, objectives.ptrs, constraints);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    result.weights = std::move(optResult.weights);
    result.iterations = optResult.iterations;
    result.finalObjective = optResult.finalObjective;
    result.converged = optResult.converged;

    Logger::info("OptimizationPipeline: " +
                 std::string(result.converged ? "converged" : "reached max iterations") +
                 " in " + std::to_string(result.iterations) + " iterations (" +
                 std::to_string(elapsed) + "s), objective=" +
                 std::to_string(result.finalObjective));

    // Forward dose computation
    auto engine = DoseEngineFactory::create("PencilBeam");
    auto dose = engine->calculateDose(dij, result.weights, doseGrid);
    result.doseResult = std::make_shared<DoseMatrix>(std::move(dose));

    // Plan analysis
    result.stats = PlanAnalysis::computeStats(
        *result.doseResult, patientData, doseGrid, config.targetDose);

    return result;
}

OptimizationPipelineResult OptimizationPipeline::run(
    const DoseInfluenceMatrix& dij,
    const OptimizationConfig& config,
    const ObjectiveProtocol& protocol,
    const PatientData& patientData,
    const Grid& doseGrid) {

    const auto& ctGrid = patientData.getGrid();
    auto objectives = ObjectiveBuilder::build(
        protocol, patientData, ctGrid, doseGrid, config.targetDose);

    return runImpl(dij, config, std::move(objectives), patientData, doseGrid);
}

OptimizationPipelineResult OptimizationPipeline::runWithObjectives(
    const DoseInfluenceMatrix& dij,
    const OptimizationConfig& config,
    BuiltObjectives objectives,
    const PatientData& patientData,
    const Grid& doseGrid) {

    return runImpl(dij, config, std::move(objectives), patientData, doseGrid);
}

} // namespace optirad
