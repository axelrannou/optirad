#pragma once

#include "DoseCalcOptions.hpp"

#include <cstddef>
#include <vector>

namespace optirad {

struct DoseObjectiveSpec {
    int type = 0;
    double threshold = 0.0, weight = 1.0, volumeFraction = 0.0;
    std::vector<size_t> indices;
};

/// Linear map from bixel weights to dose, d = A w, as the optimizers and sequencers need it.
/// Implemented by the stored sparse DoseInfluenceMatrix and by the matrix-free FluenceDoseOperator.
class IDoseOperator {
public:
    virtual ~IDoseOperator() = default;

    virtual size_t getNumVoxels() const = 0;
    virtual size_t getNumBixels() const = 0;

    virtual std::vector<double> computeDose(const std::vector<double>& weights) const = 0;

    /// grad += A^T voxelGrad
    virtual void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                            std::vector<double>& grad) const = 0;

    /// Largest dose of a single unit-weight bixel (scale for fluence bounds).
    virtual double getMaxValue() const = 0;

    /// Rows (voxels) that were computed; the others read as zero dose.
    virtual bool hasAllRows() const = 0;
    virtual bool isRowComputed(size_t voxel) const = 0;

    /// Engine settings that reproduce this operator's dose when a dose is recomputed directly.
    virtual DoseCalcOptions engineOptions() const { return {}; }
    virtual bool evaluateObjectives(const std::vector<double>&, const std::vector<DoseObjectiveSpec>&,
        double, double, double, double&, std::vector<double>&, bool = false, bool = true) const { return false; }
};

} // namespace optirad
