#include "Structure.hpp"
#include "Grid.hpp"
#include <algorithm>
#include <set>
#include <cmath>
#include <map>

namespace optirad {

// All methods are now inline in the header.
// This file can be used for future non-inline implementations
// such as complex contour processing or voxelization.

void Structure::rasterizeContours(const Grid& ctGrid) {
    m_voxelIndices.clear();
    
    if (m_contours.empty()) return;
    
    auto dims = ctGrid.getDimensions();
    auto spacing = ctGrid.getSpacing();
    auto origin = ctGrid.getOrigin();
    
    std::set<size_t> voxelSet; // Use set to avoid duplicates
    
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
    
    // Rasterize each slice
    for (const auto& [sliceIdx, contours] : contoursBySlice) {
        for (const auto* contour : contours) {
            // Use scan-line algorithm to fill polygon
            auto sliceVoxels = rasterizeContourOnSlice(*contour, ctGrid, sliceIdx);
            voxelSet.insert(sliceVoxels.begin(), sliceVoxels.end());
        }
    }
    
    m_voxelIndices.assign(voxelSet.begin(), voxelSet.end());
}

std::vector<size_t> Structure::rasterizeContourOnSlice(const Contour& contour, 
                                                        const Grid& ctGrid, 
                                                        int sliceIdx) const {
    std::vector<size_t> voxels;
    
    if (contour.points.size() < 3) return voxels;
    
    auto dims = ctGrid.getDimensions();
    auto spacing = ctGrid.getSpacing();
    auto origin = ctGrid.getOrigin();
    
    // Find bounding box of contour
    double minX = contour.points[0][0], maxX = minX;
    double minY = contour.points[0][1], maxY = minY;
    
    for (const auto& pt : contour.points) {
        minX = std::min(minX, pt[0]);
        maxX = std::max(maxX, pt[0]);
        minY = std::min(minY, pt[1]);
        maxY = std::max(maxY, pt[1]);
    }
    
    // Convert to voxel indices
    int iMin = std::max(0, static_cast<int>((minX - origin[0]) / spacing[0]));
    int iMax = std::min(static_cast<int>(dims[0]) - 1, static_cast<int>((maxX - origin[0]) / spacing[0]));
    int jMin = std::max(0, static_cast<int>((minY - origin[1]) / spacing[1]));
    int jMax = std::min(static_cast<int>(dims[1]) - 1, static_cast<int>((maxY - origin[1]) / spacing[1]));
    
    // Scan each row
    for (int j = jMin; j <= jMax; ++j) {
        double y = origin[1] + j * spacing[1];
        
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
        for (size_t i = 0; i + 1 < intersections.size(); i += 2) {
            int iStart = std::max(iMin, static_cast<int>((intersections[i] - origin[0]) / spacing[0]));
            int iEnd = std::min(iMax, static_cast<int>((intersections[i + 1] - origin[0]) / spacing[0]));
            
            for (int x = iStart; x <= iEnd; ++x) {
                size_t voxelIdx = sliceIdx * dims[0] * dims[1] + j * dims[0] + x;
                voxels.push_back(voxelIdx);
            }
        }
    }
    
    return voxels;
}

} // namespace optirad
