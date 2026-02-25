#include "SquaredOverdose.hpp"
#include <algorithm>

namespace optirad {

std::string SquaredOverdose::getName() const { return "SquaredOverdose"; }

void SquaredOverdose::setMaxDose(double dose) { m_maxDose = dose; }

double SquaredOverdose::compute(const std::vector<double>& dose) const {
    const auto& indices = getActiveIndices();
    if (indices.empty()) return 0.0;
    
    double sum = 0.0;
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        double excess = std::max(0.0, dose[idx] - m_maxDose);
        sum += excess * excess;
    }
    return m_weight * sum / indices.size();
}

std::vector<double> SquaredOverdose::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    const auto& indices = getActiveIndices();
    if (indices.empty()) return grad;
    
    double scale = 2.0 * m_weight / indices.size();
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        if (dose[idx] > m_maxDose) {
            grad[idx] = scale * (dose[idx] - m_maxDose);
        }
    }
    return grad;
}

} // namespace optirad
