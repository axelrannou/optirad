#include "SquaredOverdose.hpp"
#include <algorithm>

namespace optirad {

std::string SquaredOverdose::getName() const { return "SquaredOverdose"; }

void SquaredOverdose::setMaxDose(double dose) { m_maxDose = dose; }

double SquaredOverdose::compute(const std::vector<double>& dose) const {
    if (!m_structure) return 0.0;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) {
        return 0.0;  // No voxels to compute on
    }
    
    double sum = 0.0;
    for (size_t idx : indices) {
        double excess = std::max(0.0, dose[idx] - m_maxDose);
        sum += excess * excess;
    }
    return m_weight * sum / indices.size();
}

std::vector<double> SquaredOverdose::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    if (!m_structure) return grad;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) {
        return grad;  // No voxels to compute on, return zero gradient
    }
    
    double scale = 2.0 * m_weight / indices.size();
    for (size_t idx : indices) {
        if (dose[idx] > m_maxDose) {
            grad[idx] = scale * (dose[idx] - m_maxDose);
        }
    }
    return grad;
}

} // namespace optirad
