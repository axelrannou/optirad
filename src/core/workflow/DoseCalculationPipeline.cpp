#include "DoseCalculationPipeline.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "dose/DijSerializer.hpp"
#include "utils/Logger.hpp"

namespace optirad {

DoseCalcPipelineResult DoseCalculationPipeline::run(
    const Plan& plan,
    const Stf& stf,
    const PatientData& patientData,
    const DoseCalcPipelineOptions& options,
    DoseProgressCallback progressCallback,
    std::atomic<bool>* cancelFlag) {

    DoseCalcPipelineResult result;

    // Create dose grid
    const auto& ctGrid = patientData.getGrid();
    result.doseGrid = std::make_shared<Grid>(
        Grid::createDoseGrid(ctGrid, options.resolution));

    auto dims = result.doseGrid->getDimensions();
    Logger::info("Dose grid: " + std::to_string(dims[0]) + "x" +
                 std::to_string(dims[1]) + "x" + std::to_string(dims[2]));

    // Check cache
    if (options.useCache) {
        std::string patientName = "unknown";
        if (patientData.getPatient()) {
            patientName = patientData.getPatient()->getName();
        }
        std::string cacheFile = DijSerializer::getCacheDir() + "/" +
            DijSerializer::buildCacheKey(
                patientName,
                static_cast<int>(stf.getCount()),
                plan.getStfProperties().bixelWidth,
                options.resolution[0]);

        if (DijSerializer::exists(cacheFile)) {
            Logger::info("Loading Dij from cache: " + cacheFile);
            auto loaded = DijSerializer::load(cacheFile);
            result.dij = std::make_shared<DoseInfluenceMatrix>(std::move(loaded));
            result.cacheHit = true;
            return result;
        }
    }

    // Compute Dij
    auto engine = DoseEngineFactory::create("PencilBeam");
    if (cancelFlag) {
        engine->setCancelFlag(cancelFlag);
    }
    if (progressCallback) {
        engine->setProgressCallback(
            [&progressCallback](int current, int total, const std::string&) {
                progressCallback(current, total);
            });
    }

    DoseCalcOptions opts;
    opts.absoluteThreshold = options.absoluteThreshold;
    opts.relativeThreshold = options.relativeThreshold;
    opts.numThreads = options.numThreads;
    engine->setOptions(opts);

    auto dij = engine->calculateDij(plan, stf, patientData, *result.doseGrid);
    result.dij = std::make_shared<DoseInfluenceMatrix>(std::move(dij));

    // Save to cache
    if (options.useCache) {
        std::string patientName = "unknown";
        if (patientData.getPatient()) {
            patientName = patientData.getPatient()->getName();
        }
        std::string cacheFile = DijSerializer::getCacheDir() + "/" +
            DijSerializer::buildCacheKey(
                patientName,
                static_cast<int>(stf.getCount()),
                plan.getStfProperties().bixelWidth,
                options.resolution[0]);
        DijSerializer::save(*result.dij, cacheFile);
        Logger::info("Saved Dij to cache: " + cacheFile);
    }

    return result;
}

} // namespace optirad
