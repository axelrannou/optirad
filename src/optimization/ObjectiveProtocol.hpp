#pragma once

#include <string>
#include <vector>

namespace optirad {

/// A single objective specification in a treatment protocol.
/// Maps 1:1 to both CLI hardcoded objectives and GUI ObjectiveConfig.
struct ObjectiveSpec {
    std::string structurePattern;  // structure name or pattern
    bool exactMatch = true;        // true = exact, false = substring

    /// Objective type index matching GUI convention:
    ///   0 = SquaredDeviation, 1 = SquaredOverdose, 2 = SquaredUnderdose,
    ///   3 = MinDVH, 4 = MaxDVH
    int typeIdx = 0;

    double doseValue = 0.0;       // Gy
    double weight = 1.0;
    double volumePct = 0.0;       // % (for DVH objectives)
};

/// A treatment protocol: a named collection of objective specifications.
/// Provides a default lung-IMRT protocol used by both CLI and GUI.
struct ObjectiveProtocol {
    std::string name;
    std::vector<ObjectiveSpec> objectives;

    /// Default lung IMRT protocol matching the original hardcoded objectives.
    static ObjectiveProtocol lungIMRT(double targetDose = 66.0) {
        ObjectiveProtocol p;
        p.name = "Lung IMRT";
        p.objectives = {
            // Lungs: SquaredOverdose at MLD 20 Gy
            {"Poumon_D", true, 1, 20.0, 10.0, 0.0},
            {"Poumon_G", true, 1, 20.0, 10.0, 0.0},
            // Heart: MaxDVH V40 < 30%
            {"Coeur", true, 4, 40.0, 30.0, 30.0},
            // Esophagus: MaxDVH V50 < 30%
            {"Oesophage", true, 4, 50.0, 50.0, 30.0},
            // Spinal cord: MaxDVH V10Gy <= 0%
            {"Canal_Medullaire", true, 4, 10.0, 200.0, 0.0},
        };
        return p;
    }
};

} // namespace optirad
