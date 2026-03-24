#include "Structure.hpp"
#include "Grid.hpp"
#include <algorithm>
#include <set>
#include <cmath>
#include <map>
#ifdef OPTIRAD_HAS_TBB
#include <tbb/parallel_for.h>
#endif

namespace optirad {

// All methods are now inline in the header.
// This file can be used for future non-inline implementations
// such as complex contour processing or voxelization.

void Structure::rasterizeContours(const Grid& ctGrid) {
    // Skip rasterization if voxel indices were pre-computed (e.g., auto-generated BODY)
    if (m_preRasterized && !m_voxelIndices.empty()) {
        return;
    }
    
    m_voxelIndices.clear();
    
    if (m_contours.empty()) return;
    
    auto dims = ctGrid.getDimensions();
    auto spacing = ctGrid.getSpacing();
    auto origin = ctGrid.getOrigin();
    
    // Validate grid dimensions and spacing
    if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) {
        return;  // Invalid grid dimensions
    }
    if (std::abs(spacing[2]) < 1e-10) {
        return;  // Invalid Z spacing
    }
    
    // Group contours by Z position (slice)
    std::map<int, std::vector<const Contour*>> contoursBySlice;
    
    for (const auto& contour : m_contours) {
        if (contour.points.empty()) continue;
        
        double z = contour.zPosition;
        // Find closest slice index
        int sliceIdx = static_cast<int>(std::round((z - origin[2]) / spacing[2]));
        
        if (sliceIdx >= 0 && sliceIdx < static_cast<int>(dims[2])) {
            contoursBySlice[sliceIdx].push_back(&contour);
        }
    }
    
    if (contoursBySlice.empty()) return;

    // Prepare slice index list
    std::vector<int> sliceIndices;
    sliceIndices.reserve(contoursBySlice.size());
    for (const auto& [sliceIdx, _] : contoursBySlice) {
        sliceIndices.push_back(sliceIdx);
    }

    // Per-slice voxel buffers
    std::vector<std::vector<size_t>> sliceVoxels(sliceIndices.size());

#ifdef OPTIRAD_HAS_TBB
    // Parallel rasterization per slice
    tbb::parallel_for(size_t(0), sliceIndices.size(), [&](size_t i) {
        int sliceIdx = sliceIndices[i];
        const auto& contours = contoursBySlice.at(sliceIdx);
        std::vector<size_t> voxels;
        for (const auto* contour : contours) {
            auto v = rasterizeContourOnSlice(*contour, ctGrid, sliceIdx);
            voxels.insert(voxels.end(), v.begin(), v.end());
        }
        sliceVoxels[i] = std::move(voxels);
    });
#else
    // Sequential fallback
    for (size_t i = 0; i < sliceIndices.size(); ++i) {
        int sliceIdx = sliceIndices[i];
        const auto& contours = contoursBySlice.at(sliceIdx);
        std::vector<size_t> voxels;
        for (const auto* contour : contours) {
            auto v = rasterizeContourOnSlice(*contour, ctGrid, sliceIdx);
            voxels.insert(voxels.end(), v.begin(), v.end());
        }
        sliceVoxels[i] = std::move(voxels);
    }
#endif

    // Merge + de-duplicate
    std::vector<size_t> merged;
    size_t total = 0;
    for (const auto& v : sliceVoxels) total += v.size();
    merged.reserve(total);
    for (auto& v : sliceVoxels) {
        merged.insert(merged.end(), v.begin(), v.end());
    }

    std::sort(merged.begin(), merged.end());
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
    m_voxelIndices = std::move(merged);
}

std::vector<size_t> Structure::rasterizeContourOnSlice(const Contour& contour, 
                                                        const Grid& ctGrid, 
                                                        int sliceIdx) const {
    std::vector<size_t> voxels;
    
    if (contour.points.size() < 3) return voxels;
    
    auto dims = ctGrid.getDimensions();  // dims = [ny, nx, nz]
    auto spacing = ctGrid.getSpacing();  // spacing = [sy, sx, sz]
    auto origin = ctGrid.getOrigin();
    
    auto x_coords = ctGrid.getXCoordinates();  // Length = dims[1] = nx
    auto y_coords = ctGrid.getYCoordinates();  // Length = dims[0] = ny
    
    size_t nx = dims[1];
    size_t ny = dims[0];
    
    // Find bounding box of contour
    double minX = contour.points[0][0], maxX = minX;
    double minY = contour.points[0][1], maxY = minY;
    
    for (const auto& pt : contour.points) {
        minX = std::min(minX, pt[0]);
        maxX = std::max(maxX, pt[0]);
        minY = std::min(minY, pt[1]);
        maxY = std::max(maxY, pt[1]);
    }
    
    int jMin = 0, jMax = static_cast<int>(nx) - 1;
    for (size_t j = 0; j < nx; ++j) {
        if (x_coords[j] >= minX) { jMin = static_cast<int>(j); break; }
    }
    for (int j = static_cast<int>(nx) - 1; j >= 0; --j) {
        if (x_coords[j] <= maxX) { jMax = j; break; }
    }
    
    // Find index range for y (i indices in [0, ny-1])
    int iMin = 0, iMax = static_cast<int>(ny) - 1;
    for (size_t i = 0; i < ny; ++i) {
        if (y_coords[i] >= minY) { iMin = static_cast<int>(i); break; }
    }
    for (int i = static_cast<int>(ny) - 1; i >= 0; --i) {
        if (y_coords[i] <= maxY) { iMax = i; break; }
    }
    
    // Early exit if bounds are invalid
    if (iMin > iMax || jMin > jMax || ny == 0 || nx == 0) {
        return voxels;
    }
    
    // Scan each row in y (i index)
    for (int i = iMin; i <= iMax; ++i) {
        double y = y_coords[i];
        
        // Find intersections with contour edges at this y
        std::vector<double> intersections;
        
        for (size_t p = 0; p < contour.points.size(); ++p) {
            size_t nextP = (p + 1) % contour.points.size();
            
            double y1 = contour.points[p][1];
            double y2 = contour.points[nextP][1];
            double x1 = contour.points[p][0];
            double x2 = contour.points[nextP][0];
            
            // Check if edge crosses this scanline
            if ((y1 <= y && y < y2) || (y2 <= y && y < y1)) {
                // Compute x intersection
                double t = (y - y1) / (y2 - y1);
                double x = x1 + t * (x2 - x1);
                intersections.push_back(x);
            }
        }
        
        // Sort intersections
        std::sort(intersections.begin(), intersections.end());
        
        // Fill between pairs of intersections
        for (size_t p = 0; p + 1 < intersections.size(); p += 2) {
            // Find j indices for x coordinates
            int jStart = jMin, jEnd = jMax;
            for (int j = jMin; j <= jMax; ++j) {
                if (x_coords[j] >= intersections[p]) { jStart = j; break; }
            }
            for (int j = jMax; j >= jMin; --j) {
                if (x_coords[j] <= intersections[p + 1]) { jEnd = j; break; }
            }
            
            for (int j = jStart; j <= jEnd; ++j) {
                if (i >= 0 && i < static_cast<int>(ny) && 
                    j >= 0 && j < static_cast<int>(nx) &&
                    sliceIdx >= 0 && sliceIdx < static_cast<int>(dims[2])) {
                    size_t voxelIdx = static_cast<size_t>(i) + 
                                     static_cast<size_t>(j) * ny + 
                                     static_cast<size_t>(sliceIdx) * ny * nx;
                    voxels.push_back(voxelIdx);
                }
            }
        }
    }
    
    return voxels;
}

} // namespace optirad
