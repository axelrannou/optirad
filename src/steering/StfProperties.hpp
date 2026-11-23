#pragma once

#include <vector>
#include <array>
#include <cstddef>
#include <stdexcept>

namespace optirad {

struct StfProperties {
    std::vector<double> gantryAngles;
    std::vector<double> couchAngles;
    std::vector<double> collimatorAngles;
    size_t numOfBeams = 0;
    double bixelWidth = 7.0; // mm
    std::vector<std::array<double, 3>> isoCenters; // one per beam

    // ── Gantry angles ──

    /// Set equidistant gantry angles from start to stop with given step (degrees).
    /// If couchAngles is empty, fills it with zeros to match.
    void setGantryAngles(double start, double step, double stopExclusive) {
        gantryAngles.clear();
        for (double a = start; a < stopExclusive; a += step) {
            gantryAngles.push_back(a);
        }
        numOfBeams = gantryAngles.size();
        // Only default-fill couch angles if not already set
        if (couchAngles.empty()) {
            couchAngles.assign(numOfBeams, 0.0);
        }
        if (collimatorAngles.empty()) collimatorAngles.assign(numOfBeams, 0.0);
    }

    /// Set gantry angles from an explicit list (like matRad).
    /// If couchAngles is empty, fills it with zeros to match.
    void setGantryAngles(const std::vector<double>& angles) {
        gantryAngles = angles;
        numOfBeams = gantryAngles.size();
        if (couchAngles.empty()) {
            couchAngles.assign(numOfBeams, 0.0);
        }
        if (collimatorAngles.empty()) collimatorAngles.assign(numOfBeams, 0.0);
    }

    // ── Couch angles ──

    /// Set equidistant couch angles from start to stop with given step (degrees).
    /// When step > 0 and produces multiple couch entries, performs Cartesian product
    /// expansion: for each couch angle, ALL gantry angles are replicated (multi-arc).
    /// Result: numBeams = numGantry * numCouch.
    void setCouchAngles(double start, double step, double stopExclusive) {
        if (step <= 0.0) {
            // Single couch angle: replicate 'start' for every beam
            couchAngles.assign(numOfBeams > 0 ? numOfBeams : 1, start);
        } else {
            // Build couch range
            std::vector<double> couchRange;
            for (double a = start; a < stopExclusive; a += step) {
                couchRange.push_back(a);
            }

            if (couchRange.size() <= 1) {
                // Only one couch angle — uniform assignment
                double val = couchRange.empty() ? start : couchRange[0];
                couchAngles.assign(numOfBeams > 0 ? numOfBeams : 1, val);
            } else {
                // Cartesian product: for each couch angle, replicate all gantry angles
                std::vector<double> origGantry = gantryAngles;
                std::vector<double> origCollimator = collimatorAngles;
                gantryAngles.clear();
                couchAngles.clear();
                collimatorAngles.clear();

                for (double c : couchRange) {
                    for (size_t i = 0; i < origGantry.size(); ++i) {
                        double g = origGantry[i];
                        gantryAngles.push_back(g);
                        couchAngles.push_back(c);
                        collimatorAngles.push_back(i < origCollimator.size() ? origCollimator[i] : 0.0);
                    }
                }
                numOfBeams = gantryAngles.size();
            }
        }
    }

    /// Set couch angles from an explicit list (like matRad).
    /// Does NOT change numOfBeams or gantryAngles.
    void setCouchAngles(const std::vector<double>& angles) {
        couchAngles = angles;
    }

    /// Set collimator angles from an explicit list paired with gantry angles.
    void setCollimatorAngles(const std::vector<double>& angles) {
        collimatorAngles = angles;
    }

    /// Set a single couch angle for all beams (convenience).
    void setUniformCouchAngle(double angle) {
        couchAngles.assign(numOfBeams > 0 ? numOfBeams : 1, angle);
    }

    // ── Isocenters ──

    /// Set all isocenters to the same point
    void setUniformIsoCenter(const std::array<double, 3>& iso) {
        isoCenters.assign(numOfBeams, iso);
    }

    // ── Validation ──

    /// Ensure gantry, couch, and collimator angle lists are consistently sized.
    /// Resizes the angle lists (pad with last value or 0) if needed.
    void ensureConsistentAngles() {
        if (couchAngles.size() != gantryAngles.size()) {
            double lastCouch = couchAngles.empty() ? 0.0 : couchAngles.back();
            couchAngles.resize(gantryAngles.size(), lastCouch);
        }
        if (collimatorAngles.size() != gantryAngles.size()) {
            if (collimatorAngles.empty()) collimatorAngles.assign(gantryAngles.size(), 0.0);
            else {
                collimatorAngles.resize(gantryAngles.size(), collimatorAngles.back());
            }
        }
        numOfBeams = gantryAngles.size();
    }

    /// Returns true if angles are consistently paired.
    bool isValid() const {
        return !gantryAngles.empty() &&
               gantryAngles.size() == couchAngles.size() &&
               gantryAngles.size() == collimatorAngles.size() &&
               numOfBeams == gantryAngles.size();
    }
};

} // namespace optirad
