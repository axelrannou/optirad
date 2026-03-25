#pragma once

#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/PatientData.hpp"
#include "steering/StfProperties.hpp"
#include <memory>
#include <string>
#include <vector>
#include <array>

namespace optirad {

/// Configuration for creating a treatment plan.
struct PlanConfig {
    std::string radiationMode = "photons";
    std::string machineName = "Generic";
    int numFractions = 30;
    double bixelWidth = 7.0;

    // Gantry angles: use explicit list if non-empty, else range
    std::vector<double> gantryAngles;
    double gantryStart = 0.0;
    double gantryStep = 4.0;
    double gantryStop = 360.0;

    // Couch angles: use explicit list if non-empty, else range (step=0 means single angle)
    std::vector<double> couchAngles;
    double couchStart = 0.0;
    double couchStep = 0.0;
    double couchStop = 0.0;
};

/// Builds a treatment plan and generates STF from a PlanConfig.
/// Extracts the shared business logic from CLI createPlan() and GUI PlanningPanel.
class PlanBuilder {
public:
    struct Result {
        std::shared_ptr<Plan> plan;
        std::shared_ptr<StfProperties> stfProps;
        std::shared_ptr<Stf> stf;  // nullptr for phase-space machines
    };

    /// Build a plan and auto-generate STF for generic machines.
    /// Throws on error (machine not found, no structures, etc.).
    static Result build(const PlanConfig& config,
                        std::shared_ptr<PatientData> patientData);
};

} // namespace optirad
