#pragma once

#include "SiddonRayTracer.hpp"
#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include "core/Beam.hpp"
#include <vector>
#include <unordered_map>

namespace optirad {

/**
 * Computes radiological depth for voxels visible from a beam.
 * Ported from matRad_rayTracing.m.
 *
 * For each target voxel, traces a ray from the source through the electron density
 * volume and computes the cumulative water-equivalent path length (radiological depth).
 */
class RadDepthCalculator {
public:
    /**
     * Compute radiological depths for a set of voxels as seen from a beam.
     *
     * @param sourcePoint  Beam source point in LPS (mm)
     * @param grid         CT/dose grid
     * @param densityData  Pointer to electron density volume data
     * @param voxelIndices Linear indices of voxels to compute depths for
     * @return Map from voxel index to radiological depth (mm water-equivalent)
     */
    static std::unordered_map<size_t, double> computeRadDepths(
        const Vec3& sourcePoint,
        const Grid& grid,
        const double* densityData,
        const std::vector<size_t>& voxelIndices);

    /**
     * Compute radiological depths for all voxels along a single ray.
     * @param sourcePoint  Beam source
     * @param targetPoint  Ray target
     * @param grid         Volume grid
     * @param densityData  Electron density data
     * @return Pairs of (voxel index, radiological depth)
     */
    static std::vector<std::pair<size_t, double>> computeRayRadDepths(
        const Vec3& sourcePoint,
        const Vec3& targetPoint,
        const Grid& grid,
        const double* densityData);
};

} // namespace optirad
