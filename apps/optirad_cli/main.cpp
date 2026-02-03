#include "Logger.hpp"
#include "Config.hpp"
#include "DoseEngineFactory.hpp"
#include "OptimizerFactory.hpp"
#include "DicomImporter.hpp"
#include "Plan.hpp"
#include <iostream>

using namespace optirad;

int main(int argc, char* argv[]) {
    Logger::init();
    Logger::info("OptiRad CLI starting...");

    // Load configuration
    Config config;
    config.load("config.json");

    // Create dose engine from factory
    auto doseEngine = DoseEngineFactory::create(config.getDoseEngine());
    Logger::info("Using dose engine: " + doseEngine->getName());

    // Create optimizer from factory
    auto optimizer = OptimizerFactory::create(config.getOptimizer());
    Logger::info("Using optimizer: " + optimizer->getName());

    // Example workflow:
    // 1. Load patient data
    // DicomImporter importer;
    // importer.load("path/to/dicom");
    // auto patient = importer.getPatient();
    // auto ct = importer.getCT();
    // auto structures = importer.getStructureSet();

    // 2. Setup plan
    // Plan plan;
    // plan.addBeam(...);

    // 3. Calculate Dij
    // auto dij = doseEngine->calculateDij(plan, ct.getGrid());

    // 4. Optimize
    // auto result = optimizer->optimize(dij, objectives, constraints);

    // 5. Calculate final dose
    // auto dose = doseEngine->calculateDose(plan, ct.getGrid());

    Logger::info("OptiRad CLI finished.");
    return 0;
}
