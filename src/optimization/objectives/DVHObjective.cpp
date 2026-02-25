#include "DVHObjective.hpp"
#include <cmath>
#include <algorithm>

namespace optirad {

std::string DVHObjective::getName() const {
    return (m_type == Type::MIN_DVH) ? "MinDVH" : "MaxDVH";
}

// ─────────────────────────────────────────────────────────────────────────────
// Voxel-sorting DVH penalty (matches matRad formulation: raw sum, no 1/N)
//
// MinDVH D_v% >= d:
//   At least v% of volume should receive >= d.
//   The lowest (1-v)% are "free"; remaining voxels MUST be >= d.
//   Penalty = weight * SUM_{violating} (d - dose_i)^2
//
// MaxDVH V_d <= v%:
//   At most v% of volume should receive >= d.
//   The highest v% are "free"; remaining voxels MUST be < d.
//   Penalty = weight * SUM_{violating} (dose_i - d)^2
//
// No 1/N normalization: the raw sum scales naturally with the number of
// violating voxels, producing gradients that compete with NTO.
// ─────────────────────────────────────────────────────────────────────────────

double DVHObjective::compute(const std::vector<double>& dose) const {
    const auto& indices = getActiveIndices();
    if (indices.empty()) return 0.0;
    
    // Collect (dose, original_index) pairs and sort by dose ascending
    std::vector<std::pair<double, size_t>> sortedDoses;
    sortedDoses.reserve(indices.size());
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        sortedDoses.emplace_back(dose[idx], idx);
    }
    if (sortedDoses.empty()) return 0.0;
    
    std::sort(sortedDoses.begin(), sortedDoses.end());
    size_t N = sortedDoses.size();
    
    double obj = 0.0;
    
    if (m_type == Type::MIN_DVH) {
        size_t numFree = static_cast<size_t>(std::floor((1.0 - m_volumeFraction) * N));
        
        for (size_t i = numFree; i < N; ++i) {
            if (sortedDoses[i].first < m_doseThreshold) {
                double diff = m_doseThreshold - sortedDoses[i].first;
                obj += diff * diff;
            }
        }
    } else { // MAX_DVH
        size_t numFreeAbove = static_cast<size_t>(std::ceil(m_volumeFraction * N));
        size_t cutoff = (N > numFreeAbove) ? (N - numFreeAbove) : 0;
        
        for (size_t i = 0; i < cutoff; ++i) {
            if (sortedDoses[i].first > m_doseThreshold) {
                double diff = sortedDoses[i].first - m_doseThreshold;
                obj += diff * diff;
            }
        }
    }
    
    return m_weight * obj;
}

std::vector<double> DVHObjective::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    const auto& indices = getActiveIndices();
    if (indices.empty()) return grad;
    
    // Collect and sort by dose ascending
    std::vector<std::pair<double, size_t>> sortedDoses;
    sortedDoses.reserve(indices.size());
    for (size_t idx : indices) {
        if (idx >= dose.size()) continue;
        sortedDoses.emplace_back(dose[idx], idx);
    }
    if (sortedDoses.empty()) return grad;
    
    std::sort(sortedDoses.begin(), sortedDoses.end());
    size_t N = sortedDoses.size();
    double scale = 2.0 * m_weight;
    
    if (m_type == Type::MIN_DVH) {
        size_t numFree = static_cast<size_t>(std::floor((1.0 - m_volumeFraction) * N));
        
        for (size_t i = numFree; i < N; ++i) {
            if (sortedDoses[i].first < m_doseThreshold) {
                // gradient of (threshold - dose_i)^2 w.r.t. dose_i = -2*(threshold - dose_i)
                grad[sortedDoses[i].second] = scale * (sortedDoses[i].first - m_doseThreshold);
            }
        }
    } else { // MAX_DVH
        size_t numFreeAbove = static_cast<size_t>(std::ceil(m_volumeFraction * N));
        size_t cutoff = (N > numFreeAbove) ? (N - numFreeAbove) : 0;
        
        for (size_t i = 0; i < cutoff; ++i) {
            if (sortedDoses[i].first > m_doseThreshold) {
                // gradient of (dose_i - threshold)^2 w.r.t. dose_i = 2*(dose_i - threshold)
                grad[sortedDoses[i].second] = scale * (sortedDoses[i].first - m_doseThreshold);
            }
        }
    }
    
    return grad;
}

} // namespace optirad
