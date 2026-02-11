#include "SquaredDeviation.hpp"

namespace optirad {

std::string SquaredDeviation::getName() const { return "SquaredDeviation"; }

void SquaredDeviation::setPrescribedDose(double dose) { m_prescribedDose = dose; }

double SquaredDeviation::compute(const std::vector<double>& dose) const {
    if (!m_structure) return 0.0;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) {
        return 0.0;  // No voxels to compute on
    }
    
    double sum = 0.0;
    for (size_t idx : indices) {
        double diff = dose[idx] - m_prescribedDose;
        sum += diff * diff;
    }
    return m_weight * sum / indices.size();
}

std::vector<double> SquaredDeviation::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    if (!m_structure) return grad;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) {
        return grad;  // No voxels to compute on, return zero gradient
    }
    
    double scale = 2.0 * m_weight / indices.size();
    for (size_t idx : indices) {
        grad[idx] = scale * (dose[idx] - m_prescribedDose);
    }
    return grad;
}

} // namespace optirad
