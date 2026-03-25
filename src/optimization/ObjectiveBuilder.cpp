#include "ObjectiveBuilder.hpp"
#include "objectives/DVHObjective.hpp"
#include "objectives/SquaredOverdose.hpp"
#include "objectives/SquaredUnderdose.hpp"
#include "objectives/SquaredDeviation.hpp"
#include "geometry/StructureSet.hpp"
#include "utils/Logger.hpp"

namespace optirad {

// Helper: create an ObjectiveFunction from a spec
static std::unique_ptr<ObjectiveFunction> createObjective(const ObjectiveSpec& spec) {
    switch (spec.typeIdx) {
        case 0: {
            auto obj = std::make_unique<SquaredDeviation>();
            obj->setPrescribedDose(spec.doseValue);
            return obj;
        }
        case 1: {
            auto obj = std::make_unique<SquaredOverdose>();
            obj->setMaxDose(spec.doseValue);
            return obj;
        }
        case 2: {
            auto obj = std::make_unique<SquaredUnderdose>();
            obj->setMinDose(spec.doseValue);
            return obj;
        }
        case 3: {
            auto obj = std::make_unique<DVHObjective>();
            obj->setType(DVHObjective::Type::MIN_DVH);
            obj->setDoseThreshold(spec.doseValue);
            obj->setVolumeFraction(spec.volumePct / 100.0);
            return obj;
        }
        case 4: {
            auto obj = std::make_unique<DVHObjective>();
            obj->setType(DVHObjective::Type::MAX_DVH);
            obj->setDoseThreshold(spec.doseValue);
            obj->setVolumeFraction(spec.volumePct / 100.0);
            return obj;
        }
        default:
            return nullptr;
    }
}

// Helper: try to match a structure
static bool matchesPattern(const std::string& name, const ObjectiveSpec& spec) {
    if (spec.exactMatch) return name == spec.structurePattern;
    return name.find(spec.structurePattern) != std::string::npos;
}

BuiltObjectives ObjectiveBuilder::build(
    const ObjectiveProtocol& protocol,
    const PatientData& patientData,
    const Grid& ctGrid,
    const Grid& doseGrid,
    double targetDose,
    double targetWeight) {

    BuiltObjectives result;
    const auto* ss = patientData.getStructureSet();
    if (!ss || ss->getCount() == 0) return result;

    // Auto-add PTV objectives: MinDVH D98% for each target structure
    for (size_t i = 0; i < ss->getCount(); ++i) {
        const auto* structure = ss->getStructure(i);
        if (!structure || structure->getVoxelIndices().empty()) continue;
        if (!structure->isTarget()) continue;

        auto doseIndices = Grid::mapVoxelIndices(ctGrid, doseGrid, structure->getVoxelIndices());

        auto obj = std::make_unique<DVHObjective>();
        obj->setType(DVHObjective::Type::MIN_DVH);
        obj->setDoseThreshold(targetDose);
        obj->setVolumeFraction(0.98);
        obj->setWeight(targetWeight);
        obj->setStructure(structure);
        obj->setVoxelIndices(doseIndices);

        Logger::info("ObjectiveBuilder: Target " + structure->getName() +
                     " -> MinDVH D98% >= " + std::to_string(targetDose) +
                     " Gy (w=" + std::to_string(targetWeight) + ", " +
                     std::to_string(doseIndices.size()) + " voxels)");

        result.ptrs.push_back(obj.get());
        result.owned.push_back(std::move(obj));
    }

    // Add OAR objectives from protocol
    for (const auto& spec : protocol.objectives) {
        for (size_t i = 0; i < ss->getCount(); ++i) {
            const auto* structure = ss->getStructure(i);
            if (!structure || structure->getVoxelIndices().empty()) continue;
            if (!matchesPattern(structure->getName(), spec)) continue;

            auto doseIndices = Grid::mapVoxelIndices(ctGrid, doseGrid, structure->getVoxelIndices());
            auto obj = createObjective(spec);
            if (!obj) continue;

            obj->setWeight(spec.weight);
            obj->setStructure(structure);
            obj->setVoxelIndices(doseIndices);

            Logger::info("ObjectiveBuilder: OAR " + structure->getName() +
                         " -> type=" + std::to_string(spec.typeIdx) +
                         " dose=" + std::to_string(spec.doseValue) +
                         " Gy (w=" + std::to_string(spec.weight) + ", " +
                         std::to_string(doseIndices.size()) + " voxels)");

            result.ptrs.push_back(obj.get());
            result.owned.push_back(std::move(obj));
        }
    }

    Logger::info("ObjectiveBuilder: built " + std::to_string(result.owned.size()) + " objectives");
    return result;
}

BuiltObjectives ObjectiveBuilder::buildFromSpecs(
    const std::vector<ObjectiveSpec>& specs,
    const PatientData& patientData,
    const Grid& ctGrid,
    const Grid& doseGrid) {

    BuiltObjectives result;
    const auto* ss = patientData.getStructureSet();
    if (!ss || ss->getCount() == 0) return result;

    for (const auto& spec : specs) {
        for (size_t i = 0; i < ss->getCount(); ++i) {
            const auto* structure = ss->getStructure(i);
            if (!structure || structure->getVoxelIndices().empty()) continue;
            if (!matchesPattern(structure->getName(), spec)) continue;

            auto doseIndices = Grid::mapVoxelIndices(ctGrid, doseGrid, structure->getVoxelIndices());
            auto obj = createObjective(spec);
            if (!obj) continue;

            obj->setWeight(spec.weight);
            obj->setStructure(structure);
            obj->setVoxelIndices(doseIndices);

            result.ptrs.push_back(obj.get());
            result.owned.push_back(std::move(obj));
        }
    }

    return result;
}

} // namespace optirad
