#include "PlanBuilder.hpp"
#include "io/MachineLoader.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "utils/Logger.hpp"

namespace optirad {

PlanBuilder::Result PlanBuilder::build(const PlanConfig& config,
                                       std::shared_ptr<PatientData> patientData) {
    Result result;

    auto plan = std::make_shared<Plan>();
    plan->setName("TreatmentPlan");
    plan->setRadiationMode(config.radiationMode);
    plan->setNumOfFractions(config.numFractions);
    plan->setPatientData(patientData);

    // Load machine
    plan->setMachine(MachineLoader::load(config.radiationMode, config.machineName));

    // Configure STF properties
    StfProperties stfProps;
    if (!config.gantryAngles.empty()) {
        stfProps.setGantryAngles(config.gantryAngles);
    } else {
        stfProps.setGantryAngles(config.gantryStart, config.gantryStep, config.gantryStop);
    }

    if (!config.couchAngles.empty()) {
        stfProps.setCouchAngles(config.couchAngles);
    } else if (config.couchStep > 0.0) {
        stfProps.setCouchAngles(config.couchStart, config.couchStep, config.couchStop);
    } else {
        stfProps.setUniformCouchAngle(config.couchStart);
    }

    stfProps.ensureConsistentAngles();
    stfProps.bixelWidth = config.bixelWidth;

    // Compute isocenter from target structures
    auto iso = plan->computeIsoCenter();
    stfProps.setUniformIsoCenter(iso);

    plan->setStfProperties(stfProps);
    result.plan = plan;

    // Auto-generate STF for generic (non-phase-space) machines
    if (!plan->getMachine().isPhaseSpace()) {
        const auto& stfP = plan->getStfProperties();
        std::array<double, 3> isoRef = {0.0, 0.0, 0.0};
        if (!stfP.isoCenters.empty()) {
            isoRef = stfP.isoCenters[0];
        }

        PhotonIMRTStfGenerator generator(0.0, 360.0, 360.0, stfP.bixelWidth, isoRef);
        generator.setMachine(plan->getMachine());
        generator.setRadiationMode(plan->getRadiationMode());
        generator.setGantryAngles(stfP.gantryAngles);
        generator.setCouchAngles(stfP.couchAngles);

        if (patientData && patientData->hasValidCT() && patientData->hasStructures()) {
            const auto* ct = patientData->getCTVolume();
            const auto* structureSet = patientData->getStructureSet();
            const auto& grid = ct->getGrid();

            generator.setGrid(grid);
            generator.setStructureSet(*structureSet);
            generator.setCTResolution(grid.getSpacing());
        }

        result.stfProps = generator.generate();
        result.stf = std::make_shared<Stf>(generator.generateStf());
    }

    return result;
}

} // namespace optirad
