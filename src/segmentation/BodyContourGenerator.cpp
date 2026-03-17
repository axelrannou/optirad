#include "BodyContourGenerator.hpp"
#include "utils/Logger.hpp"
#include <queue>
#include <algorithm>
#include <cmath>
#include <mutex>

#ifdef OPTIRAD_HAS_TBB
#include <tbb/parallel_for.h>
#endif

namespace optirad {

std::unique_ptr<Structure> BodyContourGenerator::generate(const Volume<int16_t>& ctVolume,
                                                           int16_t huThreshold) {
    const Grid& grid = ctVolume.getGrid();
    auto dims = grid.getDimensions();
    // dims[0] = nRows (Y direction), dims[1] = nCols (X direction), dims[2] = nSlices
    size_t nRows   = dims[0];
    size_t nCols   = dims[1];
    size_t nSlices = dims[2];

    if (nRows == 0 || nCols == 0 || nSlices == 0) {
        Logger::error("BodyContourGenerator: invalid CT dimensions");
        return nullptr;
    }

    Logger::info("Auto-generating BODY contour from CT (HU threshold: " + 
                 std::to_string(huThreshold) + ")");

    const int16_t* ctData = ctVolume.data();

    // Per-slice results
    std::vector<Contour> contours(nSlices);
    std::vector<std::vector<size_t>> sliceVoxelIndices(nSlices);
    std::vector<bool> sliceHasBody(nSlices, false);

#ifdef OPTIRAD_HAS_TBB
    Logger::info("Using TBB parallel body contour generation");
    tbb::parallel_for(size_t(0), nSlices, [&](size_t k) {
        auto mask = createSliceMask(ctData, nRows, nCols, k, huThreshold);
        
        // Check if slice has any body voxels
        bool hasBody = std::any_of(mask.begin(), mask.end(), [](bool v) { return v; });
        if (!hasBody) return;

        // Collect voxel indices directly from the mask (no contour roundtrip)
        sliceVoxelIndices[k] = collectSliceVoxelIndices(mask, nRows, nCols, k);
        
        // Extract contour polygon for GUI visualization
        contours[k] = extractSliceContour(mask, nRows, nCols, grid, k);
        sliceHasBody[k] = true;
    });
#else
    Logger::info("Using sequential body contour generation");
    for (size_t k = 0; k < nSlices; ++k) {
        auto mask = createSliceMask(ctData, nRows, nCols, k, huThreshold);
        
        bool hasBody = std::any_of(mask.begin(), mask.end(), [](bool v) { return v; });
        if (!hasBody) continue;

        sliceVoxelIndices[k] = collectSliceVoxelIndices(mask, nRows, nCols, k);
        contours[k] = extractSliceContour(mask, nRows, nCols, grid, k);
        sliceHasBody[k] = true;
    }
#endif

    // Assemble the BODY structure
    auto bodyStructure = std::make_unique<Structure>();
    bodyStructure->setName("BODY");
    bodyStructure->setType("EXTERNAL");
    bodyStructure->setROINumber(0);
    bodyStructure->setColor(139, 90, 43);  // Brown
    bodyStructure->setPriority(5);

    // Merge all voxel indices
    size_t totalVoxels = 0;
    size_t contourCount = 0;
    for (size_t k = 0; k < nSlices; ++k) {
        totalVoxels += sliceVoxelIndices[k].size();
        if (sliceHasBody[k]) ++contourCount;
    }

    if (totalVoxels == 0) {
        Logger::warn("BodyContourGenerator: no body voxels found");
        return nullptr;
    }

    // Set voxel indices directly (bypasses contour→rasterize roundtrip)
    std::vector<size_t> allVoxelIndices;
    allVoxelIndices.reserve(totalVoxels);
    for (size_t k = 0; k < nSlices; ++k) {
        allVoxelIndices.insert(allVoxelIndices.end(), 
                               sliceVoxelIndices[k].begin(), 
                               sliceVoxelIndices[k].end());
    }
    bodyStructure->setVoxelIndices(allVoxelIndices);
    bodyStructure->setPreRasterized(true);

    // Add contours for GUI visualization
    for (size_t k = 0; k < nSlices; ++k) {
        if (sliceHasBody[k] && !contours[k].points.empty()) {
            bodyStructure->addContour(contours[k]);
        }
    }

    Logger::info("  - BODY (EXTERNAL, " + std::to_string(contourCount) + 
                 " contours, " + std::to_string(totalVoxels) + " voxels, RGB: 139/90/43)");

    return bodyStructure;
}

std::vector<bool> BodyContourGenerator::createSliceMask(const int16_t* ctData,
                                                         size_t nRows, size_t nCols,
                                                         size_t sliceIdx,
                                                         int16_t huThreshold) {
    size_t sliceSize = nRows * nCols;
    size_t sliceOffset = sliceIdx * sliceSize;
    std::vector<bool> mask(sliceSize, false);

    // Threshold — column-major layout: data[row + col * nRows + slice * nRows * nCols]
    for (size_t col = 0; col < nCols; ++col) {
        for (size_t row = 0; row < nRows; ++row) {
            size_t idx = row + col * nRows;
            if (ctData[sliceOffset + idx] > huThreshold) {
                mask[idx] = true;
            }
        }
    }

    // Flood-fill exterior to create solid body envelope
    floodFillExterior(mask, nRows, nCols);

    // Keep only the largest connected component (patient body)
    // This removes couch, table, and small artifacts
    keepLargestComponent(mask, nRows, nCols);

    return mask;
}

void BodyContourGenerator::floodFillExterior(std::vector<bool>& mask, size_t nRows, size_t nCols) {
    // Mark exterior air pixels reachable from borders.
    // After flood-fill, any pixel that is NOT body AND NOT exterior = internal cavity → set to true.
    
    size_t total = nRows * nCols;
    std::vector<bool> exterior(total, false);
    std::queue<std::array<size_t, 2>> queue;

    // Helper: mask index for (row, col) in column-major layout
    auto idx = [nRows](size_t row, size_t col) -> size_t {
        return row + col * nRows;
    };

    // Seed from border pixels that are NOT body
    // Top and bottom rows
    for (size_t col = 0; col < nCols; ++col) {
        size_t topIdx = idx(0, col);
        size_t botIdx = idx(nRows - 1, col);
        if (!mask[topIdx] && !exterior[topIdx]) {
            exterior[topIdx] = true;
            queue.push({0, col});
        }
        if (!mask[botIdx] && !exterior[botIdx]) {
            exterior[botIdx] = true;
            queue.push({nRows - 1, col});
        }
    }
    // Left and right columns
    for (size_t row = 0; row < nRows; ++row) {
        size_t leftIdx = idx(row, 0);
        size_t rightIdx = idx(row, nCols - 1);
        if (!mask[leftIdx] && !exterior[leftIdx]) {
            exterior[leftIdx] = true;
            queue.push({row, 0});
        }
        if (!mask[rightIdx] && !exterior[rightIdx]) {
            exterior[rightIdx] = true;
            queue.push({row, nCols - 1});
        }
    }

    // BFS flood-fill (4-connected)
    const int dr[] = {-1, 1, 0, 0};
    const int dc[] = {0, 0, -1, 1};

    while (!queue.empty()) {
        auto [cr, cc] = queue.front();
        queue.pop();

        for (int d = 0; d < 4; ++d) {
            int nr = static_cast<int>(cr) + dr[d];
            int nc = static_cast<int>(cc) + dc[d];

            if (nr < 0 || nr >= static_cast<int>(nRows) || 
                nc < 0 || nc >= static_cast<int>(nCols)) continue;

            size_t nIdx = idx(static_cast<size_t>(nr), static_cast<size_t>(nc));
            if (!mask[nIdx] && !exterior[nIdx]) {
                exterior[nIdx] = true;
                queue.push({static_cast<size_t>(nr), static_cast<size_t>(nc)});
            }
        }
    }

    // Fill internal cavities: any pixel that is NOT body AND NOT exterior → make it body
    for (size_t i = 0; i < total; ++i) {
        if (!mask[i] && !exterior[i]) {
            mask[i] = true;
        }
    }
}

void BodyContourGenerator::keepLargestComponent(std::vector<bool>& mask, size_t nRows, size_t nCols) {
    size_t total = nRows * nCols;
    // Label array: 0 = unlabeled/background
    std::vector<int> labels(total, 0);
    int currentLabel = 0;
    int largestLabel = 0;
    size_t largestSize = 0;

    auto idx = [nRows](size_t row, size_t col) -> size_t {
        return row + col * nRows;
    };

    const int dr[] = {-1, 1, 0, 0};
    const int dc[] = {0, 0, -1, 1};

    // BFS-based connected component labeling
    for (size_t col = 0; col < nCols; ++col) {
        for (size_t row = 0; row < nRows; ++row) {
            size_t i = idx(row, col);
            if (!mask[i] || labels[i] != 0) continue;

            // New component found — BFS
            ++currentLabel;
            size_t componentSize = 0;
            std::queue<std::array<size_t, 2>> q;
            labels[i] = currentLabel;
            q.push({row, col});

            while (!q.empty()) {
                auto [cr, cc] = q.front();
                q.pop();
                ++componentSize;

                for (int d = 0; d < 4; ++d) {
                    int nr = static_cast<int>(cr) + dr[d];
                    int nc = static_cast<int>(cc) + dc[d];
                    if (nr < 0 || nr >= static_cast<int>(nRows) ||
                        nc < 0 || nc >= static_cast<int>(nCols)) continue;

                    size_t ni = idx(static_cast<size_t>(nr), static_cast<size_t>(nc));
                    if (mask[ni] && labels[ni] == 0) {
                        labels[ni] = currentLabel;
                        q.push({static_cast<size_t>(nr), static_cast<size_t>(nc)});
                    }
                }
            }

            if (componentSize > largestSize) {
                largestSize = componentSize;
                largestLabel = currentLabel;
            }
        }
    }

    // Zero out everything except the largest component
    if (largestLabel > 0) {
        for (size_t i = 0; i < total; ++i) {
            if (mask[i] && labels[i] != largestLabel) {
                mask[i] = false;
            }
        }
    }
}

std::vector<size_t> BodyContourGenerator::collectSliceVoxelIndices(const std::vector<bool>& mask,
                                                                    size_t nRows, size_t nCols,
                                                                    size_t sliceIdx) {
    std::vector<size_t> indices;
    indices.reserve(mask.size() / 4);  // Rough estimate

    size_t sliceStride = nRows * nCols;
    for (size_t col = 0; col < nCols; ++col) {
        for (size_t row = 0; row < nRows; ++row) {
            size_t maskIdx = row + col * nRows;
            if (mask[maskIdx]) {
                // The mask uses Volume layout (transposed from DICOM row-major).
                // Volume::at(row, col) = DICOM(DICOM_row=col, DICOM_col=row)
                // Flat index must use DICOM order: i + j * nRows where i=DICOM_row, j=DICOM_col
                // Therefore: voxelIdx = col + row * nRows + slice * nRows * nCols + 1
                size_t voxelIdx = col + row * nRows + sliceIdx * sliceStride + 1;
                indices.push_back(voxelIdx);
            }
        }
    }

    return indices;
}

Contour BodyContourGenerator::extractSliceContour(const std::vector<bool>& mask,
                                                   size_t nRows, size_t nCols,
                                                   const Grid& grid,
                                                   size_t sliceIdx) {
    // Row-scanning approach: for each row, find the leftmost and rightmost body column.
    // Build a closed polygon from left edges (top→bottom) + right edges (bottom→top).
    // This produces a clean, non-self-intersecting polygon for visualization.

    auto maskIdx = [nRows](size_t row, size_t col) -> size_t {
        return row + col * nRows;
    };

    // Find row range with body pixels
    size_t minRow = nRows, maxRow = 0;
    for (size_t row = 0; row < nRows; ++row) {
        for (size_t col = 0; col < nCols; ++col) {
            if (mask[maskIdx(row, col)]) {
                minRow = std::min(minRow, row);
                maxRow = std::max(maxRow, row);
                break;  // Only need to find if this row has any body pixel
            }
        }
    }

    if (minRow > maxRow) return {};  // No body pixels

    // For each row in the body range, find left and right body edges
    struct RowEdge { size_t row; size_t leftCol; size_t rightCol; };
    std::vector<RowEdge> edges;
    edges.reserve(maxRow - minRow + 1);

    for (size_t row = minRow; row <= maxRow; ++row) {
        size_t leftCol = nCols, rightCol = 0;
        for (size_t col = 0; col < nCols; ++col) {
            if (mask[maskIdx(row, col)]) {
                if (col < leftCol) leftCol = col;
                if (col > rightCol) rightCol = col;
            }
        }
        if (leftCol <= rightCol) {
            edges.push_back({row, leftCol, rightCol});
        }
    }

    if (edges.empty()) return {};

    // Build polygon: left edges top→bottom, then right edges bottom→top
    // Points placed at pixel edges (±0.5 offset) for proper enclosure
    std::vector<std::array<double, 2>> polyPoints;
    polyPoints.reserve(edges.size() * 2 + 2);

    // Left edge (top to bottom) — offset col by -0.5
    for (const auto& e : edges) {
        polyPoints.push_back({static_cast<double>(e.row), static_cast<double>(e.leftCol) - 0.5});
    }
    // Right edge (bottom to top) — offset col by +0.5
    for (auto it = edges.rbegin(); it != edges.rend(); ++it) {
        polyPoints.push_back({static_cast<double>(it->row), static_cast<double>(it->rightCol) + 0.5});
    }

    // Subsample if polygon has too many points (keep up to ~600 for lightweight display)
    if (polyPoints.size() > 600) {
        size_t step = polyPoints.size() / 300;
        if (step < 2) step = 2;
        std::vector<std::array<double, 2>> subsampled;
        subsampled.reserve(polyPoints.size() / step + 1);
        for (size_t i = 0; i < polyPoints.size(); i += step) {
            subsampled.push_back(polyPoints[i]);
        }
        polyPoints = std::move(subsampled);
    }

    // Convert to patient LPS coordinates
    Contour contour;
    contour.points.reserve(polyPoints.size());

    for (const auto& [pRow, pCol] : polyPoints) {
        // The mask uses flat index row + col * nRows which equals Volume::at(row, col, slice).
        // Due to DICOM row-major → Volume column-major transposition, Volume index (row, col)
        // corresponds to DICOM pixel (DICOM_row=col, DICOM_col=row).
        // voxelToPatient expects (DICOM_row, DICOM_col, slice) = (col, row, slice).
        Vec3 ijk = {pCol, pRow, static_cast<double>(sliceIdx)};
        Vec3 lps = grid.voxelToPatient(ijk);
        contour.points.push_back(lps);
    }

    if (!contour.points.empty()) {
        contour.zPosition = contour.points[0][2];
    }

    return contour;
}

} // namespace optirad
