#pragma once

#include "../ObjectiveFunction.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Grid.hpp"
#include <memory>

namespace optirad {

/// Eclipse-style Normal Tissue Objective. Doses are in % of the prescription; distances in mm.
struct NtoSettings {
    bool enabled = true;
    double priority = 100.0;
    double distanceFromTargetMm = 1.0;    // the limit stays at the start dose up to this distance
    double startDosePct = 105.0;
    double endDosePct = 60.0;
    double fallOffPerMm = 0.05;
};

/// Upper limit on the dose of normal tissue (body minus all targets) that falls off with the distance x from the
/// target border, as in Eclipse:
///   limit(x) = start                                                  for x <= distance
///   limit(x) = start * e^(-k (x - distance)) + end * (1 - e^(-k (x - distance)))   beyond
/// Cost = priority * mean over the normal tissue voxels of max(0, dose - limit)^2: like every objective it is
/// normalized by its voxel count, so the priority does not depend on the body size or the dose grid.
class NormalTissueObjective : public ObjectiveFunction {
public:
    NormalTissueObjective(std::vector<size_t> voxels, std::vector<double> limits, std::string summary = "");

    /// Builds the objective on the dose grid, or returns null when it is disabled or there is no target.
    /// With `dose`, voxels without a dose row (body voxels skipped by the Dij) are left out.
    static std::unique_ptr<NormalTissueObjective> build(const NtoSettings& settings, double prescriptionDose,
                                                        const StructureSet* structures, const Grid& ctGrid,
                                                        const Grid& doseGrid, const IDoseOperator* dose = nullptr);

    /// Dose limit (Gy) at `distanceMm` from the target border.
    static double limit(const NtoSettings& settings, double prescriptionDose, double distanceMm);

    std::string getName() const override { return "NTO"; }
    std::string describe() const override;
    bool gpuDescription(DoseObjectiveSpec& output) const override {
        output = {1, 0.0, m_weight, 0.0, getActiveIndices(), m_limits};
        return true;
    }
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

    const std::vector<double>& limits() const { return m_limits; }
    const std::vector<size_t>& voxels() const { return getActiveIndices(); }

private:
    std::vector<double> m_limits;   // per voxel, parallel to the mapped indices
    std::string m_summary;
};

} // namespace optirad
