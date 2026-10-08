#pragma once

#include "Plan.hpp"
#include "Grid.hpp"
#include "Stf.hpp"
#include "PatientData.hpp"
#include "DoseMatrix.hpp"
#include "DoseInfluenceMatrix.hpp"
#include "DoseCalcOptions.hpp"
#include <string>
#include <functional>
#include <atomic>
#include <stdexcept>

namespace optirad {

/// Progress callback: (currentBeam, totalBeams, message)
using DoseCalcProgressCallback = std::function<void(int, int, const std::string&)>;

class IDoseEngine {
public:
    virtual ~IDoseEngine() = default;

    virtual std::string getName() const = 0;

    /// Calculate dose influence matrix (Dij) for optimization.
    /// @param plan        Treatment plan (machine data, fractions, etc.)
    /// @param stf         Steering file with beam/ray geometry
    /// @param patientData Patient CT + electron density + structures
    /// @param doseGrid    Grid for dose computation (may differ from CT grid)
    /// @return Sparse Dij matrix [numDoseVoxels x totalNumBixels]
    virtual DoseInfluenceMatrix calculateDij(
        const Plan& plan,
        const Stf& stf,
        const PatientData& patientData,
        const Grid& doseGrid) = 0;

    /// Calculate forward dose from Dij and weights.
    /// @param dij     Pre-computed dose influence matrix
    /// @param weights Beamlet fluence weights
    /// @param grid    Dose grid for result
    /// @return DoseMatrix with computed dose
    virtual DoseMatrix calculateDose(
        const DoseInfluenceMatrix& dij,
        const std::vector<double>& weights,
        const Grid& grid) = 0;

    /// Compute the dose of a weight vector directly (no Dij kept in memory), over all structure voxels.
    /// Needed for the full dose cube when the Dij skips external voxels, and for deliverable dose.
    virtual DoseMatrix calculateDoseDirect(
        const Plan& plan,
        const Stf& stf,
        const PatientData& patientData,
        const Grid& doseGrid,
        const std::vector<double>& weights) {
        (void)plan; (void)stf; (void)patientData; (void)doseGrid; (void)weights;
        throw std::runtime_error("calculateDoseDirect is not supported by engine " + getName());
    }

    /// Set progress callback for GUI/CLI progress reporting
    void setProgressCallback(DoseCalcProgressCallback cb) { m_progressCallback = std::move(cb); }

    /// Set cancellation flag (atomic, thread-safe)
    void setCancelFlag(std::atomic<bool>* flag) { m_cancelFlag = flag; }

    /// Set options for dose calculation (thresholds, parallelism).
    void setOptions(const DoseCalcOptions& opts) { m_options = opts; }

    /// Get current options.
    const DoseCalcOptions& getOptions() const { return m_options; }

protected:
    DoseCalcProgressCallback m_progressCallback;
    std::atomic<bool>* m_cancelFlag = nullptr;
    DoseCalcOptions m_options;
};

} // namespace optirad
