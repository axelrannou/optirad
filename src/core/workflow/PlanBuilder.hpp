#pragma once

#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/PatientData.hpp"
#include "core/ArcPlan.hpp"
#include "steering/StfProperties.hpp"
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <optional>

namespace optirad {

/// Configuration for creating a treatment plan.
struct PlanConfig {
    std::string radiationMode = "photons";
    std::string machineName = "TrueBeamSTx_6X";
    int numFractions = 30;
    double bixelWidth = 3.0;

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
    std::vector<double> collimatorAngles;

    // Delivery technique. For VMAT the gantry/couch fields above are ignored: the fluence-optimization
    // beams are the sector centres of the arcs in `vmat`.
    DeliveryMode deliveryMode = DeliveryMode::VMAT;
    VmatSettings vmat;

    // Isocenter (patient LPS mm). When unset: the imported RTPLAN isocenter if the patient has one,
    // otherwise the centre of the target structures.
    std::optional<std::array<double, 3>> isocenter;

    // Dose-grid resolution (mm). When unset: the imported RTDOSE spacing if any, otherwise the Plan default.
    std::optional<std::array<double, 3>> doseGridResolution;

    // Whole-course prescription (Gy). When unset the Plan keeps its default.
    std::optional<double> prescribedDose;
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

    /// Isocenter precedence: `config.isocenter`, then the imported RTPLAN isocenter, then
    /// `plan.computeIsoCenter()` (target centre), so planning works with or without an imported plan.
    static std::array<double, 3> resolveIsocenter(const PlanConfig& config, const Plan& plan,
                                                  const PatientData* patientData);
};

} // namespace optirad
