#include "RadDepthCalculator.hpp"
#include <cmath>
#include <algorithm>

namespace optirad {

std::unordered_map<size_t, double> RadDepthCalculator::computeRadDepths(
    const Vec3& sourcePoint,
    const Grid& grid,
    const double* densityData,
    const std::vector<size_t>& voxelIndices)
{
    std::unordered_map<size_t, double> radDepths;
    radDepths.reserve(voxelIndices.size());

    auto dims = grid.getDimensions();
    size_t nx = dims[0], ny = dims[1];

    // For each target voxel, trace a ray from source through that voxel
    for (size_t idx : voxelIndices) {
        // Convert linear index to ijk
        size_t iz = idx / (nx * ny);
        size_t rem = idx % (nx * ny);
        size_t iy = rem / nx;
        size_t ix = rem % nx;

        // Convert to LPS coordinates
        Vec3 ijk = {static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)};
        Vec3 targetLPS = grid.voxelToPatient(ijk);

        // Trace ray from source to target voxel
        auto rayDepths = SiddonRayTracer::traceRadDepth(sourcePoint, targetLPS, grid, densityData);

        // Find the depth entry for this specific voxel (or closest to it)
        // The last entry in the trace is closest to the target
        if (!rayDepths.empty()) {
            // Look for exact voxel match
            bool found = false;
            for (const auto& [vIdx, depth] : rayDepths) {
                if (vIdx == idx) {
                    radDepths[idx] = depth;
                    found = true;
                    break;
                }
            }
            // If not found (rounding), use the last entry
            if (!found) {
                radDepths[idx] = rayDepths.back().second;
            }
        }
    }

    return radDepths;
}

std::vector<std::pair<size_t, double>> RadDepthCalculator::computeRayRadDepths(
    const Vec3& sourcePoint,
    const Vec3& targetPoint,
    const Grid& grid,
    const double* densityData)
{
    return SiddonRayTracer::traceRadDepth(sourcePoint, targetPoint, grid, densityData);
}

} // namespace optirad
