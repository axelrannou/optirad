#include "SquaredUnderdose.hpp"
#include <algorithm>

namespace optirad {

std::string SquaredUnderdose::getName() const { return "SquaredUnderdose"; }

void SquaredUnderdose::setMinDose(double dose) { m_minDose = dose; }

double SquaredUnderdose::compute(const std::vector<double>& dose) const {
    const auto& indices = getActiveIndices();
    if (indices.empty()) return 0.0;
    
    double sum = 0.0;
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        double deficit = std::max(0.0, m_minDose - dose[idx]);
        sum += deficit * deficit;
    }
    return m_weight * sum / indices.size();
}

std::vector<double> SquaredUnderdose::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    const auto& indices = getActiveIndices();
    if (indices.empty()) return grad;
    
    double scale = 2.0 * m_weight / indices.size();
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        if (dose[idx] < m_minDose) {
            grad[idx] = -scale * (m_minDose - dose[idx]);
        }
    }
    return grad;
}

} // namespace optirad
