#pragma once

#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include <vector>
#include <cstddef>

namespace optirad {

/**
 * Result of a single ray trace through a volume grid.
 */
struct RayTraceResult {
    std::vector<size_t> voxelIndices;          ///< Linear indices of traversed voxels
    std::vector<double> intersectionLengths;   ///< Path length through each voxel (mm)
    std::vector<double> densities;             ///< Density (electron density) at each voxel
    double totalDistance = 0.0;                ///< Total geometric distance from source to target (mm)
};

/**
 * Siddon ray tracer (Siddon 1985) for computing ray-voxel intersections
 * through a 3D volume grid. Used for radiological depth computation and SSD calculation.
 *
 * The algorithm computes parametric intersections of a ray (from source to target)
 * with all grid planes, then determines which voxels are traversed and the
 * intersection lengths within each voxel.
 */
class SiddonRayTracer {
public:
    /**
     * Trace a ray through a volume grid.
     * @param source       Source point in patient (LPS) coordinates (mm)
     * @param target       Target point in patient (LPS) coordinates (mm)
     * @param grid         The volume grid (defines dimensions, spacing, origin)
     * @param densityData  Pointer to density volume data (electron density), can be nullptr
     * @return RayTraceResult with traversed voxels, lengths, and densities
     */
    static RayTraceResult trace(
        const Vec3& source,
        const Vec3& target,
        const Grid& grid,
        const double* densityData = nullptr);

    /**
     * Trace a ray and compute cumulative radiological depth per traversed voxel.
     * Uses midpoint rule: radDepth[i] = sum(l[0..i-1] * rho[0..i-1]) + l[i]*rho[i]/2
     * @param source       Source point in LPS (mm)
     * @param target       Target point in LPS (mm)
     * @param grid         Volume grid
     * @param densityData  Pointer to electron density data
     * @return Vector of (voxelIndex, cumulativeRadDepth) pairs
     */
    static std::vector<std::pair<size_t, double>> traceRadDepth(
        const Vec3& source,
        const Vec3& target,
        const Grid& grid,
        const double* densityData);
};

} // namespace optirad
