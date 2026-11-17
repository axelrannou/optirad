#pragma once

#include "ObjectiveProtocol.hpp"
#include "ObjectiveFunction.hpp"
#include "objectives/NormalTissueObjective.hpp"
#include "core/PatientData.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <vector>
#include <string>

namespace optirad {

/// Built objectives ready for the optimizer.
struct BuiltObjectives {
    std::vector<std::unique_ptr<ObjectiveFunction>> owned;
    std::vector<ObjectiveFunction*> ptrs;  // raw pointers for IOptimizer API
};

/// Builds ObjectiveFunction instances from an ObjectiveProtocol + patient structures.
/// Extracts the shared logic from CLI optimize() and GUI OptimizationPanel thread.
class ObjectiveBuilder {
public:
    /// Build objectives from a protocol, adding PTV targets from the structure set.
    /// Target structures are auto-detected and get MinDVH D98% >= Rx, MaxDVH 0% <= 105% Rx and
    /// SquaredDeviation(Rx) objectives (the GUI "Recommended" template).
    /// @param protocol     OAR objectives specification
    /// @param patientData  Patient with structure set
    /// @param ctGrid       CT grid (for voxel index mapping)
    /// @param doseGrid     Dose grid (for voxel index mapping)
    /// @param targetDose   Prescribed dose for PTV objectives (Gy); required, throws if not > 0
    /// @param targetWeight Weight for PTV objectives
    static BuiltObjectives build(
        const ObjectiveProtocol& protocol,
        const PatientData& patientData,
        const Grid& ctGrid,
        const Grid& doseGrid,
        double targetDose,
        double targetWeight = 100.0);

    /// Build objectives from a flat list of ObjectiveSpec (without auto-PTV detection).
    /// Used by the GUI when the user has manually specified all objectives.
    static BuiltObjectives buildFromSpecs(
        const std::vector<ObjectiveSpec>& specs,
        const PatientData& patientData,
        const Grid& ctGrid,
        const Grid& doseGrid);

    /// Appends the Eclipse-style normal tissue objective (no-op when disabled or without a target).
    /// Returns true when an NTO was added.
    static bool addNormalTissueObjective(
        BuiltObjectives& objectives,
        const NtoSettings& settings,
        double prescriptionDose,
        const PatientData& patientData,
        const Grid& doseGrid,
        const IDoseOperator* dose = nullptr);
};

} // namespace optirad
