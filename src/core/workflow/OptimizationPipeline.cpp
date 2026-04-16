#include "OptimizationPipeline.hpp"
#include "optimization/OptimizerFactory.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/Constraint.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "utils/Logger.hpp"
#include <chrono>
#include <unordered_map>
#include <cmath>

namespace optirad {

/// Build a 2D neighbor map from STF beam geometry.
/// For each bixel (global index), stores indices of its 4-connected neighbors
/// in the same beam's BEV grid at ±bixelWidth in X and Z.
static std::vector<std::vector<int>> buildNeighborMap(const Stf& stf) {
    size_t totalBixels = stf.getTotalNumOfBixels();
    std::vector<std::vector<int>> neighbors(totalBixels);

    size_t globalOffset = 0;
    for (size_t bi = 0; bi < stf.getCount(); ++bi) {
        const auto* beam = stf.getBeam(bi);
        if (!beam || beam->getNumOfRays() == 0) {
            globalOffset += beam ? beam->getTotalNumOfBixels() : 0;
            continue;
        }

        double bw = beam->getBixelWidth();
        size_t nRays = beam->getNumOfRays();

        // Build a map from quantized BEV (col, row) to global bixel index
        struct PairHash {
            size_t operator()(const std::pair<int,int>& p) const {
                return std::hash<long long>()(static_cast<long long>(p.first) * 100000LL + p.second);
            }
        };
        std::unordered_map<std::pair<int,int>, size_t, PairHash> posToIdx;

        // First pass: find grid origin
        double xMin = 1e30, zMin = 1e30;
        for (size_t ri = 0; ri < nRays; ++ri) {
            const auto& pos = beam->getRay(ri)->getRayPosBev();
            if (pos[0] < xMin) xMin = pos[0];
            if (pos[2] < zMin) zMin = pos[2];
        }

        // Second pass: map positions to grid indices
        for (size_t ri = 0; ri < nRays; ++ri) {
            const auto& pos = beam->getRay(ri)->getRayPosBev();
            int col = static_cast<int>(std::round((pos[0] - xMin) / bw));
            int row = static_cast<int>(std::round((pos[2] - zMin) / bw));
            posToIdx[{col, row}] = globalOffset + ri;
        }

        // Third pass: link neighbors
        for (size_t ri = 0; ri < nRays; ++ri) {
            const auto& pos = beam->getRay(ri)->getRayPosBev();
            int col = static_cast<int>(std::round((pos[0] - xMin) / bw));
            int row = static_cast<int>(std::round((pos[2] - zMin) / bw));
            size_t gIdx = globalOffset + ri;

            static const int dx[] = {1, -1, 0, 0};
            static const int dz[] = {0, 0, 1, -1};
            for (int d = 0; d < 4; ++d) {
                auto it = posToIdx.find({col + dx[d], row + dz[d]});
                if (it != posToIdx.end()) {
                    neighbors[gIdx].push_back(static_cast<int>(it->second));
                }
            }
        }

        globalOffset += beam->getTotalNumOfBixels();
    }

    return neighbors;
}

static OptimizationPipelineResult runImpl(
    const DoseInfluenceMatrix& dij,
    const OptimizationConfig& config,
    BuiltObjectives objectives,
    const PatientData& patientData,
    const Grid& doseGrid,
    IterationCallback iterCallback,
    const Stf* stf) {

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

        // Spatial smoothing
        if (config.spatialSmoothingWeight > 0.0 && stf) {
            auto neighbors = buildNeighborMap(*stf);
            lbfgs->setSpatialSmoothing(config.spatialSmoothingWeight, neighbors);
        }

        // Regularization
        if (config.l2RegWeight > 0.0) lbfgs->setL2Regularization(config.l2RegWeight);
        if (config.l1RegWeight > 0.0) lbfgs->setL1Regularization(config.l1RegWeight);
    }

    if (iterCallback) {
        optimizer->setIterationCallback(std::move(iterCallback));
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
    const Grid& doseGrid,
    const Stf* stf) {

    const auto& ctGrid = patientData.getGrid();
    auto objectives = ObjectiveBuilder::build(
        protocol, patientData, ctGrid, doseGrid, config.targetDose);

    return runImpl(dij, config, std::move(objectives), patientData, doseGrid, nullptr, stf);
}

OptimizationPipelineResult OptimizationPipeline::runWithObjectives(
    const DoseInfluenceMatrix& dij,
    const OptimizationConfig& config,
    BuiltObjectives objectives,
    const PatientData& patientData,
    const Grid& doseGrid,
    IterationCallback iterCallback,
    const Stf* stf) {

    return runImpl(dij, config, std::move(objectives), patientData, doseGrid, std::move(iterCallback), stf);
}

} // namespace optirad
