#include "SiddonRayTracer.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace optirad {

RayTraceResult SiddonRayTracer::trace(
    const Vec3& source,
    const Vec3& target,
    const Grid& grid,
    const double* densityData)
{
    RayTraceResult result;

    auto dims = grid.getDimensions();
    auto spacing = grid.getSpacing();
    auto origin = grid.getOrigin();

    size_t nx = dims[0], ny = dims[1], nz = dims[2];

    // Ray direction vector
    Vec3 dir = vecSub(target, source);
    double d12 = norm(dir);
    result.totalDistance = d12;

    if (d12 < 1e-10) return result;

    // Grid plane positions (boundaries of voxels)
    // Plane i is at origin[d] + (i - 0.5) * spacing[d] for i = 0..N
    // Or equivalently: first plane at origin - spacing/2, last at origin + (N-0.5)*spacing

    // Compute parametric intersections alpha_min and alpha_max for the volume
    // The volume spans from origin - spacing/2 to origin + (N-0.5)*spacing in each axis
    double alphaMin = 0.0;
    double alphaMax = 1.0;

    // Plane coordinates: planes at x = origin[0] + (i - 0.5)*spacing[0] for i = 0..nx
    // In parametric form: source + alpha * dir crosses plane at alpha = (plane - source) / dir

    // Arrays to hold parametric values for each axis
    std::vector<double> alphas;
    alphas.reserve(nx + ny + nz + 6);

    for (int axis = 0; axis < 3; ++axis) {
        size_t N = dims[axis];
        double sp = spacing[axis];
        double orig = origin[axis];
        double d = dir[axis];
        double s = source[axis];

        if (std::abs(d) < 1e-12) continue; // Ray parallel to this axis

        // First and last planes
        double planeFirst = orig - sp * 0.5;
        double planeLast = orig + (static_cast<double>(N) - 0.5) * sp;

        double a1 = (planeFirst - s) / d;
        double aN = (planeLast - s) / d;

        double aLow = std::min(a1, aN);
        double aHigh = std::max(a1, aN);

        alphaMin = std::max(alphaMin, aLow);
        alphaMax = std::min(alphaMax, aHigh);

        if (alphaMin >= alphaMax) {
            // Ray misses the volume
            return result;
        }

        // Compute all plane crossings between alphaMin and alphaMax
        double aStart = std::min(a1, aN);
        double step = sp / std::abs(d);

        // Find the first plane crossing >= alphaMin
        int iStart, iEnd;
        if (d > 0) {
            iStart = static_cast<int>(std::ceil((s + alphaMin * d - orig + sp * 0.5) / sp));
            iEnd   = static_cast<int>(std::floor((s + alphaMax * d - orig + sp * 0.5) / sp));
        } else {
            iStart = static_cast<int>(std::ceil((s + alphaMax * d - orig + sp * 0.5) / sp));
            iEnd   = static_cast<int>(std::floor((s + alphaMin * d - orig + sp * 0.5) / sp));
        }

        iStart = std::max(iStart, 0);
        iEnd   = std::min(iEnd, static_cast<int>(N));

        for (int i = iStart; i <= iEnd; ++i) {
            double plane = orig + (static_cast<double>(i) - 0.5) * sp;
            double a = (plane - s) / d;
            if (a >= alphaMin && a <= alphaMax) {
                alphas.push_back(a);
            }
        }
    }

    if (alphas.empty()) return result;

    // Add alphaMin and alphaMax
    alphas.push_back(alphaMin);
    alphas.push_back(alphaMax);

    // Sort and remove duplicates
    std::sort(alphas.begin(), alphas.end());
    alphas.erase(std::unique(alphas.begin(), alphas.end(),
        [](double a, double b) { return std::abs(a - b) < 1e-12; }),
        alphas.end());

    // Compute voxel crossings
    for (size_t i = 0; i + 1 < alphas.size(); ++i) {
        double aMid = 0.5 * (alphas[i] + alphas[i + 1]);
        double len = (alphas[i + 1] - alphas[i]) * d12;

        if (len < 1e-12) continue;

        // Midpoint position
        Vec3 mid = {
            source[0] + aMid * dir[0],
            source[1] + aMid * dir[1],
            source[2] + aMid * dir[2]
        };

        // Convert to voxel indices
        Vec3 ijk = grid.patientToVoxel(mid);
        int ix = static_cast<int>(std::round(ijk[0]));
        int iy = static_cast<int>(std::round(ijk[1]));
        int iz = static_cast<int>(std::round(ijk[2]));

        // Bounds check
        if (ix < 0 || ix >= static_cast<int>(nx) ||
            iy < 0 || iy >= static_cast<int>(ny) ||
            iz < 0 || iz >= static_cast<int>(nz)) {
            continue;
        }

        // Linear index (row-major to match DICOM pixel storage: row*cols + col + slice*rows*cols)
        // ix = row index (dims[0] = Rows), iy = col index (dims[1] = Columns)
        size_t linearIdx = static_cast<size_t>(ix) * ny +
                           static_cast<size_t>(iy) +
                           static_cast<size_t>(iz) * nx * ny;

        result.voxelIndices.push_back(linearIdx);
        result.intersectionLengths.push_back(len);

        if (densityData) {
            result.densities.push_back(densityData[linearIdx]);
        } else {
            result.densities.push_back(1.0); // Water-equivalent default
        }
    }

    return result;
}

std::vector<std::pair<size_t, double>> SiddonRayTracer::traceRadDepth(
    const Vec3& source,
    const Vec3& target,
    const Grid& grid,
    const double* densityData)
{
    auto traceResult = trace(source, target, grid, densityData);

    std::vector<std::pair<size_t, double>> radDepths;
    radDepths.reserve(traceResult.voxelIndices.size());

    double cumRadDepth = 0.0;

    for (size_t i = 0; i < traceResult.voxelIndices.size(); ++i) {
        double l = traceResult.intersectionLengths[i];
        double rho = traceResult.densities[i];

        // Midpoint rule: add half of current voxel first
        double midDepth = cumRadDepth + 0.5 * l * rho;
        radDepths.emplace_back(traceResult.voxelIndices[i], midDepth);

        // Accumulate full crossing
        cumRadDepth += l * rho;
    }

    return radDepths;
}

} // namespace optirad
