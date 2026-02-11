#pragma once

#include <vector>
#include <array>
#include <cstddef>

namespace optirad {

struct StfProperties {
    std::vector<double> gantryAngles;
    std::vector<double> couchAngles;
    size_t numOfBeams = 0;
    double bixelWidth = 7.0; // mm
    std::vector<std::array<double, 3>> isoCenters; // one per beam

    /// Helper: set equidistant gantry angles from start to stop with given step (degrees)
    void setGantryAngles(double start, double step, double stopExclusive) {
        gantryAngles.clear();
        for (double a = start; a < stopExclusive; a += step) {
            gantryAngles.push_back(a);
        }
        numOfBeams = gantryAngles.size();
        couchAngles.assign(numOfBeams, 0.0);
    }

    /// Set all isocenters to the same point
    void setUniformIsoCenter(const std::array<double, 3>& iso) {
        isoCenters.assign(numOfBeams, iso);
    }
};

} // namespace optirad
