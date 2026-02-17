#pragma once

#include <vector>
#include <array>
#include <set>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <mutex>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

/// Performs 3D morphological dilation on a set of voxel indices.
///
/// Expands the target voxel set by the given margin (in mm) using
/// 26-connectivity. The expansion is constrained to stay within the
/// patient surface (union of all structure voxels).
///
/// @param targetVoxelIndices   1-based linear indices of target voxels
/// @param allVoxelIndices      1-based linear indices of all structure voxels (patient surface)
/// @param dimensions           Grid dimensions {ny, nx, nz}
/// @param spacing              Grid spacing {dy, dx, dz} in mm
/// @param margin               Margin to expand in mm {mx, my, mz}
/// @return                     Expanded set of 1-based linear voxel indices
inline std::vector<size_t> dilateVoxels(
    const std::vector<size_t>& targetVoxelIndices,
    const std::set<size_t>& allVoxelIndices,
    const std::array<size_t, 3>& dimensions,
    const std::array<double, 3>& spacing,
    const std::array<double, 3>& margin)
{
    size_t ny = dimensions[0];
    size_t nx = dimensions[1];
    size_t nz = dimensions[2];

    // Number of voxels to expand in each axis
    int voxelMarginY = static_cast<int>(std::round(margin[0] / spacing[0]));
    int voxelMarginX = static_cast<int>(std::round(margin[1] / spacing[1]));
    int voxelMarginZ = static_cast<int>(std::round(margin[2] / spacing[2]));
    int maxIter = std::max({voxelMarginY, voxelMarginX, voxelMarginZ});

    if (maxIter <= 0) {
        return targetVoxelIndices;
    }

#ifdef _OPENMP
    std::cout << "  Performing 3D dilation with OpenMP (" << omp_get_max_threads() 
              << " threads, " << maxIter << " iterations)..." << std::flush;
#else
    std::cout << "  Performing 3D dilation sequentially (" << maxIter << " iterations)..." << std::flush;
#endif

    // Build the enlarged set starting with the target
    std::set<size_t> enlarged(targetVoxelIndices.begin(), targetVoxelIndices.end());

    for (int cnt = 1; cnt <= maxIter; ++cnt) {
        // The frontier: voxels added in the previous iteration (or initial set for cnt=1)
        // On first iteration, frontier is the full set; on subsequent iterations,
        // it's the difference from the previous pass
        std::vector<size_t> frontier(enlarged.begin(), enlarged.end());

        // Decompose frontier voxels to subscripts, exclude border voxels
        struct Sub { size_t row; size_t col; size_t slice; };
        std::vector<Sub> subs;
        subs.reserve(frontier.size());

        for (size_t idx1 : frontier) {
            size_t idx0 = idx1 - 1;
            size_t row   = idx0 % ny;
            size_t col   = (idx0 / ny) % nx;
            size_t slice = idx0 / (ny * nx);

            // Exclude border voxels (can't expand beyond grid)
            if (row == 0 || row == ny - 1 ||
                col == 0 || col == nx - 1 ||
                slice == 0 || slice == nz - 1) {
                continue;
            }

            subs.push_back({row, col, slice});
        }

        // Which axes are active this iteration
        int dy = (voxelMarginY >= cnt) ? 1 : 0;
        int dx = (voxelMarginX >= cnt) ? 1 : 0;
        int dz = (voxelMarginZ >= cnt) ? 1 : 0;

        // Collect new voxels in thread-local buffers, then merge
        std::vector<std::vector<size_t>> threadBuffers;
#ifdef _OPENMP
        int nThreads = omp_get_max_threads();
        threadBuffers.resize(nThreads);

        // 26-connectivity expansion in parallel
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& localBuffer = threadBuffers[tid];
            
            #pragma omp for collapse(3) schedule(dynamic)
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    for (int dk = -1; dk <= 1; ++dk) {
                        if (di == 0 && dj == 0 && dk == 0) continue;

                        for (const auto& s : subs) {
                            size_t newRow   = s.row   + di * dy;
                            size_t newCol   = s.col   + dj * dx;
                            size_t newSlice = s.slice + dk * dz;

                            // Convert back to 1-based linear index
                            size_t newIdx1 = newRow + newCol * ny + newSlice * ny * nx + 1;

                            // Constrain to patient surface
                            if (allVoxelIndices.count(newIdx1)) {
                                localBuffer.push_back(newIdx1);
                            }
                        }
                    }
                }
            }
        }
        
        // Merge thread-local buffers into enlarged set
        for (const auto& buffer : threadBuffers) {
            enlarged.insert(buffer.begin(), buffer.end());
        }
#else
        // Sequential fallback: 26-connectivity expansion
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    if (di == 0 && dj == 0 && dk == 0) continue;

                    for (const auto& s : subs) {
                        size_t newRow   = s.row   + di * dy;
                        size_t newCol   = s.col   + dj * dx;
                        size_t newSlice = s.slice + dk * dz;

                        // Convert back to 1-based linear index (column-major: row + col*ny + slice*ny*nx)
                        size_t newIdx1 = newRow + newCol * ny + newSlice * ny * nx + 1;

                        // Constrain to patient surface
                        if (allVoxelIndices.count(newIdx1)) {
                            enlarged.insert(newIdx1);
                        }
    std::cout << " done\n" << std::flush;
                    }
                }
            }
        }
#endif
    }

    return std::vector<size_t>(enlarged.begin(), enlarged.end());
}

} // namespace optirad
