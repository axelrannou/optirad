#pragma once

#include "Structure.hpp"
#include <vector>
#include <string>

namespace optirad {

class ObjectiveFunction {
public:
    virtual ~ObjectiveFunction() = default;

    virtual std::string getName() const = 0;

    // Compute objective value given dose distribution
    virtual double compute(const std::vector<double>& dose) const = 0;

    // Compute gradient with respect to dose
    virtual std::vector<double> gradient(const std::vector<double>& dose) const = 0;

    void setWeight(double weight);
    double getWeight() const;

    void setStructure(const Structure* structure);

    /// Set pre-mapped voxel indices (e.g., dose-grid indices).
    /// When set, these are used instead of m_structure->getVoxelIndices().
    void setVoxelIndices(const std::vector<size_t>& indices);

protected:
    /// Return the voxel indices to iterate over.
    /// Uses m_mappedIndices if non-empty, otherwise falls back to m_structure.
    const std::vector<size_t>& getActiveIndices() const;

    double m_weight = 1.0;
    const Structure* m_structure = nullptr;
    std::vector<size_t> m_mappedIndices;  ///< dose-grid indices (set by caller)
};

} // namespace optirad
