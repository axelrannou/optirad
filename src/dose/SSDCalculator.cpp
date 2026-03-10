#include "SSDCalculator.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace optirad {

double SSDCalculator::computeSSD(
    const Vec3& source,
    const Vec3& target,
    const Grid& grid,
    const double* densityData,
    double densityThreshold)
{
    auto traceResult = SiddonRayTracer::trace(source, target, grid, densityData);

    if (traceResult.voxelIndices.empty()) {
        return -1.0; // Ray misses volume
    }

    // Find first voxel above density threshold (surface entry point)
    double cumulativeLength = 0.0;
    for (size_t i = 0; i < traceResult.voxelIndices.size(); ++i) {
        if (traceResult.densities[i] > densityThreshold) {
            // SSD = distance from source to this voxel entry
            // Approximate as cumulative length up to entry of this voxel
            return cumulativeLength + norm(vecSub(source, source)); // distance from source to volume entry + cumulative
        }
        cumulativeLength += traceResult.intersectionLengths[i];
    }

    // Didn't find surface - compute based on first crossing alpha
    // The ray entered the volume at some distance from source
    Vec3 dir = vecSub(target, source);
    double d12 = norm(dir);

    // First voxel entry is approximately at the start of the trace
    // Better: compute actual entry distance
    auto dims = grid.getDimensions();
    auto spacing = grid.getSpacing();
    auto origin = grid.getOrigin();

    // Compute alphaMin (entry into volume)
    double alphaMin = 0.0;
    double alphaMax = 1.0;

    for (int axis = 0; axis < 3; ++axis) {
        double sp = spacing[axis];
        double orig = origin[axis];
        double d = dir[axis];
        double s = source[axis];
        size_t N = dims[axis];

        if (std::abs(d) < 1e-12) continue;

        double planeFirst = orig - sp * 0.5;
        double planeLast = orig + (static_cast<double>(N) - 0.5) * sp;

        double a1 = (planeFirst - s) / d;
        double aN = (planeLast - s) / d;

        alphaMin = std::max(alphaMin, std::min(a1, aN));
        alphaMax = std::min(alphaMax, std::max(a1, aN));
    }

    // Walk through trace to find surface
    double distFromEntry = alphaMin * d12;
    for (size_t i = 0; i < traceResult.densities.size(); ++i) {
        if (traceResult.densities[i] > densityThreshold) {
            return distFromEntry;
        }
        distFromEntry += traceResult.intersectionLengths[i];
    }

    return -1.0; // No surface found
}

std::vector<double> SSDCalculator::computeBeamSSDs(
    const Vec3& sourcePoint,
    const std::vector<Vec3>& rayTargets,
    const Grid& grid,
    const double* densityData,
    double densityThreshold)
{
    size_t numRays = rayTargets.size();
    std::vector<double> ssds(numRays, -1.0);

    // Compute SSD for each ray
    for (size_t i = 0; i < numRays; ++i) {
        ssds[i] = computeSSD(sourcePoint, rayTargets[i], grid, densityData, densityThreshold);
    }

    // Fix missing SSDs from nearest valid ray (matRad: matRad_closestNeighbourSSD)
    // Find rays with valid SSDs
    std::vector<size_t> validIndices;
    for (size_t i = 0; i < numRays; ++i) {
        if (ssds[i] > 0.0) {
            validIndices.push_back(i);
        }
    }

    if (validIndices.empty()) {
        // Fallback: use SAD as default SSD  
        std::fill(ssds.begin(), ssds.end(), 1000.0);
        return ssds;
    }

    // For missing SSDs, find nearest valid ray by target position distance
    for (size_t i = 0; i < numRays; ++i) {
        if (ssds[i] > 0.0) continue;

        double minDist = std::numeric_limits<double>::max();
        size_t nearestIdx = validIndices[0];

        for (size_t vi : validIndices) {
            double dist = norm(vecSub(rayTargets[i], rayTargets[vi]));
            if (dist < minDist) {
                minDist = dist;
                nearestIdx = vi;
            }
        }

        ssds[i] = ssds[nearestIdx];
    }

    return ssds;
}

} // namespace optirad
