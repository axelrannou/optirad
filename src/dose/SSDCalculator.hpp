#pragma once

#include "SiddonRayTracer.hpp"
#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include <vector>

namespace optirad {

/**
 * Computes Source-to-Surface Distance (SSD) for each ray in a beam.
 * Ported from matRad_computeSSD.m.
 *
 * SSD is the distance from the source point along the ray direction
 * to the first voxel where electron density exceeds a threshold (skin surface).
 */
class SSDCalculator {
public:
    /**
     * Compute SSD for a single ray.
     * @param source       Source point in LPS coordinates (mm)
     * @param target       Target point in LPS coordinates (mm)
     * @param grid         Volume grid
     * @param densityData  Pointer to electron density data
     * @param densityThreshold  Minimum density to be considered "surface" (default 0.05)
     * @return SSD in mm, or -1 if ray misses the volume
     */
    static double computeSSD(
        const Vec3& source,
        const Vec3& target,
        const Grid& grid,
        const double* densityData,
        double densityThreshold = 0.05);

    /**
     * Compute SSD for all rays in a beam.
     * Rays that miss the volume get their SSD from the nearest valid ray.
     * @param sourcePoint  Beam source point in LPS
     * @param rayTargets   Target points for each ray in LPS
     * @param grid         Volume grid
     * @param densityData  Pointer to electron density data
     * @param densityThreshold  Surface threshold
     * @return Vector of SSDs (one per ray)
     */
    static std::vector<double> computeBeamSSDs(
        const Vec3& sourcePoint,
        const std::vector<Vec3>& rayTargets,
        const Grid& grid,
        const double* densityData,
        double densityThreshold = 0.05);
};

} // namespace optirad
