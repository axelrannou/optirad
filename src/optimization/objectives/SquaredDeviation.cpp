#include "SquaredDeviation.hpp"

namespace optirad {

std::string SquaredDeviation::getName() const { return "SquaredDeviation"; }

void SquaredDeviation::setPrescribedDose(double dose) { m_prescribedDose = dose; }

double SquaredDeviation::compute(const std::vector<double>& dose) const {
    const auto& indices = getActiveIndices();
    if (indices.empty()) return 0.0;
    
    double sum = 0.0;
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        double diff = dose[idx] - m_prescribedDose;
        sum += diff * diff;
    }
    return m_weight * sum / indices.size();
}

std::vector<double> SquaredDeviation::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    const auto& indices = getActiveIndices();
    if (indices.empty()) return grad;
    
    double scale = 2.0 * m_weight / indices.size();
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        grad[idx] = scale * (dose[idx] - m_prescribedDose);
    }
    return grad;
}

} // namespace optirad
