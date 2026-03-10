#pragma once

#include "DoseMatrix.hpp"
#include "core/PatientData.hpp"
#include "geometry/Grid.hpp"
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <iostream>

namespace optirad {

/// Per-structure dose statistics for plan analysis.
struct StructureDoseStats {
    std::string name;
    std::string type;        // e.g., "TARGET", "OAR", "EXTERNAL"
    size_t numVoxels = 0;    // number of dose-grid voxels

    double minDose = 0.0;
    double maxDose = 0.0;
    double meanDose = 0.0;
    double stdDose = 0.0;

    // Dx% values (dose covering x% of volume)
    double d2 = 0.0;     // D2%  — near-max dose
    double d5 = 0.0;     // D5%
    double d50 = 0.0;    // D50% — median dose
    double d95 = 0.0;    // D95%
    double d98 = 0.0;    // D98% — near-min dose

    // VxGy values (% volume receiving >= x Gy)
    double v20 = 0.0;    // V20Gy [%]
    double v40 = 0.0;    // V40Gy [%]
    double v50 = 0.0;    // V50Gy [%]
    double v60 = 0.0;    // V60Gy [%]

    // Quality indices (only for targets)
    double conformityIndex = 0.0;  // CI = V_ref_in_target / V_target
    double homogeneityIndex = 0.0; // HI = (D2 - D98) / D50

    /// Sorted dose values for DVH curve computation.
    std::vector<double> sortedDoses;
};

/// Compute DVH curve points from sorted dose values.
/// Returns (dose, volume%) pairs suitable for plotting.
struct DVHCurveData {
    std::string structureName;
    std::vector<float> doses;     // x-axis: dose [Gy]
    std::vector<float> volumes;   // y-axis: volume [%]
};

/**
 * Plan quality analysis utility — computes per-structure statistics and DVH curves.
 * Equivalent to matRad's matRad_planAnalysis function.
 */
class PlanAnalysis {
public:
    /// Compute all statistics for all structures.
    /// @param dose  The dose result (dose-grid sized)
    /// @param patient  Patient data with structure set
    /// @param doseGrid  The dose calculation grid
    /// @param prescribedDose  Prescribed dose for target (used for CI calculation)
    static std::vector<StructureDoseStats> computeStats(
        const DoseMatrix& dose,
        const PatientData& patient,
        const Grid& doseGrid,
        double prescribedDose = 60.0);

    /// Compute DVH curve data for each structure.
    /// @param stats  Pre-computed statistics (must have sortedDoses populated)
    /// @param maxDose  Maximum dose for x-axis range
    /// @param numBins  Number of dose bins for the DVH curve
    static std::vector<DVHCurveData> computeDVHCurves(
        const std::vector<StructureDoseStats>& stats,
        double maxDose,
        int numBins = 200);

    /// Print plan analysis to stream (matRad-style output).
    static void print(const std::vector<StructureDoseStats>& stats,
                      std::ostream& os = std::cout);

    /// Compute Dx% — dose covering x% of the volume.
    /// sortedDoses must be sorted ascending.
    static double computeDx(const std::vector<double>& sortedDoses, double percentile);

    /// Compute VxGy — percentage of volume receiving >= doseThreshold.
    static double computeVx(const std::vector<double>& sortedDoses, double doseThreshold);
};

} // namespace optirad
