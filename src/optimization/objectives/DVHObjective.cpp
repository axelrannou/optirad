#include "DVHObjective.hpp"
#include <cmath>
#include <algorithm>

namespace optirad {

std::string DVHObjective::getName() const { return "DVHObjective"; }

double DVHObjective::compute(const std::vector<double>& dose) const {
    if (!m_structure) return 0.0;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) return 0.0;
    
    // Smooth Heaviside approximation for differentiability
    // H(x) ≈ 0.5 * (1 + tanh(k * x)) where k controls sharpness
    constexpr double k = 10.0;
    
    double sum = 0.0;
    for (size_t idx : indices) {
        double diff = dose[idx] - m_doseThreshold;
        double h = 0.5 * (1.0 + std::tanh(k * diff));  // ~1 if dose > threshold
        sum += h;
    }
    
    double actualFraction = sum / indices.size();
    double violation = 0.0;
    
    if (m_type == Type::MIN_DOSE) {
        // Penalize if actual fraction < target fraction
        violation = std::max(0.0, m_volumeFraction - actualFraction);
    } else {
        // Penalize if actual fraction > target fraction  
        violation = std::max(0.0, actualFraction - m_volumeFraction);
    }
    
    return m_weight * violation * violation;
}

std::vector<double> DVHObjective::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    if (!m_structure) return grad;
    
    auto indices = m_structure->getVoxelIndices();
    if (indices.empty()) return grad;
    
    constexpr double k = 10.0;
    double n = static_cast<double>(indices.size());
    
    // Compute actual volume fraction and its derivative
    double sum = 0.0;
    for (size_t idx : indices) {
        double diff = dose[idx] - m_doseThreshold;
        sum += 0.5 * (1.0 + std::tanh(k * diff));
    }
    double actualFraction = sum / n;
    
    double violation = 0.0;
    double sign = 0.0;
    
    if (m_type == Type::MIN_DOSE) {
        violation = m_volumeFraction - actualFraction;
        if (violation <= 0) return grad;
        sign = -1.0;  // Increasing dose decreases violation
    } else {
        violation = actualFraction - m_volumeFraction;
        if (violation <= 0) return grad;
        sign = 1.0;
    }
    
    // d/d(dose_i) of violation^2 = 2 * violation * d(violation)/d(dose_i)
    // d(actualFraction)/d(dose_i) = (1/n) * d(H_i)/d(dose_i)
    // d(H)/d(x) = 0.5 * k * sech^2(k*x)
    
    double scale = 2.0 * m_weight * violation * sign / n;
    
    for (size_t idx : indices) {
        double diff = dose[idx] - m_doseThreshold;
        double sech = 1.0 / std::cosh(k * diff);
        double dH = 0.5 * k * sech * sech;
        grad[idx] = scale * dH;
    }
    
    return grad;
}

} // namespace optirad
