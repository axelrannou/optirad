#include "Machine.hpp"

#include <algorithm>
#include <cmath>

namespace optirad {

double MachineData::relativeFluence(double radiusMm) const {
    if (!useCustomPrimaryFluence || primaryFluence.empty()) return 1.0;
    if (radiusMm <= primaryFluence.front()[0]) return primaryFluence.front()[1];
    if (radiusMm >= primaryFluence.back()[0]) return primaryFluence.back()[1];
    const auto upper = std::lower_bound(primaryFluence.begin(), primaryFluence.end(), radiusMm,
        [](const std::array<double, 2>& entry, double r) { return entry[0] < r; });
    const auto& hi = *upper;
    const auto& lo = *(upper - 1);
    const double t = hi[0] > lo[0] ? (radiusMm - lo[0]) / (hi[0] - lo[0]) : 0.0;
    return lo[1] + t * (hi[1] - lo[1]);
}

double MachineData::depthDoseCorrection(double radDepthMm) const {
    if (extraAttenuation == 0.0) return 1.0;
    return std::exp(-extraAttenuation * (radDepthMm - extraAttenuationRefDepth));
}

} // namespace optirad
