#pragma once

#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/PatientData.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseCalcOptions.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <functional>
#include <atomic>
#include <array>
#include <string>

namespace optirad {

/// Options for the dose calculation pipeline.
struct DoseCalcPipelineOptions {
    std::array<double, 3> resolution = {2.5, 2.5, 2.5};
    bool useCache = true;
    double absoluteThreshold = 1e-6;
    double relativeThreshold = 0.01;   // fraction (0.01 = 1%)
    int numThreads = 0;                // 0 = all
};

/// Result of the dose calculation pipeline.
struct DoseCalcPipelineResult {
    std::shared_ptr<DoseInfluenceMatrix> dij;
    std::shared_ptr<Grid> doseGrid;
    bool cacheHit = false;
};

/// Progress callback: (currentBeam, totalBeams)
using DoseProgressCallback = std::function<void(int, int)>;

/// Extracts the shared Dij computation pipeline from CLI doseCalc() and GUI PlanningPanel.
/// Handles: grid creation → cache lookup → compute → cache save.
class DoseCalculationPipeline {
public:
    /// Run the full pipeline. Throws on error.
    /// cancelFlag can be nullptr if cancellation is not needed.
    static DoseCalcPipelineResult run(
        const Plan& plan,
        const Stf& stf,
        const PatientData& patientData,
        const DoseCalcPipelineOptions& options,
        DoseProgressCallback progressCallback = nullptr,
        std::atomic<bool>* cancelFlag = nullptr);
};

} // namespace optirad
