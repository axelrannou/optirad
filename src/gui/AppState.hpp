#pragma once

#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "steering/StfProperties.hpp"
#include <memory>

namespace optirad {

/// Shared application state for the GUI pipeline.
/// Mirrors the CLI AppState: load DICOM → create Plan → generate STF.
struct GuiAppState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;
    std::shared_ptr<Stf> stf;

    // Workflow state
    bool dicomLoaded() const { return patientData != nullptr; }
    bool planCreated() const { return plan != nullptr; }
    bool stfGenerated() const { return stf != nullptr && !stf->isEmpty(); }

    // Reset downstream state when upstream changes
    void resetPlan() { plan.reset(); stfProps.reset(); stf.reset(); }
    void resetStf() { stfProps.reset(); stf.reset(); }
};

} // namespace optirad
