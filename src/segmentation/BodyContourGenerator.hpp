#pragma once

#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <vector>
#include <cstdint>

namespace optirad {

/**
 * Automatic BODY/EXTERNAL contour generator from CT volume.
 *
 * When no BODY (EXTERNAL) structure is present in the RT-STRUCT,
 * this class generates one by:
 *   1. Thresholding the CT at a configurable HU value (default: -300)
 *   2. Flood-filling from slice borders to identify exterior air
 *   3. Filling internal cavities (lungs, airways) to produce a solid body envelope
 *   4. Connected-component analysis: keep only the largest component per slice
 *      (removes couch, table, and small artifacts)
 *   5. Computing voxel indices directly from the mask (for dose calculation)
 *   6. Extracting per-slice contour polygons from the mask (for GUI visualization)
 *
 * The resulting Structure has type "EXTERNAL", name "BODY", pre-computed
 * voxel indices, and contour polygons for display.
 */
class BodyContourGenerator {
public:
    /// Default HU threshold for body segmentation
    static constexpr int16_t DEFAULT_HU_THRESHOLD = -300;

    /**
     * Generate a BODY structure from a CT volume.
     *
     * @param ctVolume  The CT volume in Hounsfield Units
     * @param huThreshold  HU threshold (voxels > threshold are body). Default: -400
     * @return A Structure with type "EXTERNAL", pre-computed voxel indices, and contour polygons
     */
    static std::unique_ptr<Structure> generate(const Volume<int16_t>& ctVolume,
                                                int16_t huThreshold = DEFAULT_HU_THRESHOLD);

private:
    /**
     * Create a binary body mask for a single axial slice.
     * Thresholds at huThreshold, then flood-fills exterior from borders
     * to produce a solid body envelope (internal cavities filled).
     *
     * @param ctData       Raw CT data pointer
     * @param nRows, nCols Slice dimensions (dims[0], dims[1])
     * @param sliceIdx     Slice index (k)
     * @param huThreshold  HU threshold
     * @return 2D binary mask (nRows * nCols), indexed as [col * nRows + row], true = body
     */
    static std::vector<bool> createSliceMask(const int16_t* ctData,
                                              size_t nRows, size_t nCols,
                                              size_t sliceIdx,
                                              int16_t huThreshold);

    /**
     * Flood-fill exterior air from slice borders.
     * Marks all non-body pixels reachable from borders as exterior.
     * Remaining non-body pixels are internal cavities → filled.
     *
     * @param mask       In/out binary mask (modified in place)
     * @param nRows, nCols Slice dimensions
     */
    static void floodFillExterior(std::vector<bool>& mask, size_t nRows, size_t nCols);

    /**
     * Keep only the largest 4-connected component in the binary mask.
     * Performs connected-component labeling via BFS, then zeroes out
     * all components except the largest one (the patient body).
     * This removes couch, table, and small disconnected artifacts.
     *
     * @param mask       In/out binary mask (modified in place)
     * @param nRows, nCols Slice dimensions
     */
    static void keepLargestComponent(std::vector<bool>& mask, size_t nRows, size_t nCols);

    /**
     * Collect 1-based voxel indices from the mask for a single slice.
     * Uses the same column-major indexing as Structure::rasterizeContourOnSlice:
     *   flatIndex = row + col * nRows + sliceIdx * nRows * nCols + 1
     *
     * @param mask       Binary mask for the slice
     * @param nRows, nCols Slice dimensions
     * @param sliceIdx   Slice index
     * @return Vector of 1-based voxel indices
     */
    static std::vector<size_t> collectSliceVoxelIndices(const std::vector<bool>& mask,
                                                         size_t nRows, size_t nCols,
                                                         size_t sliceIdx);

    /**
     * Extract a contour polygon from the mask for GUI visualization.
     * Uses row-scanning to find left/right body edges per row,
     * producing a clean, non-self-intersecting polygon.
     *
     * @param mask       Binary mask for the slice
     * @param nRows, nCols Slice dimensions
     * @param grid       CT grid for coordinate conversion
     * @param sliceIdx   Slice index
     * @return Contour with patient-space points and zPosition set
     */
    static Contour extractSliceContour(const std::vector<bool>& mask,
                                        size_t nRows, size_t nCols,
                                        const Grid& grid,
                                        size_t sliceIdx);
};

} // namespace optirad
