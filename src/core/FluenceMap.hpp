#pragma once

#include "Beam.hpp"
#include "Machine.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace optirad {

/// 2D fluence map extracted from optimizer weights for a single beam.
/// Grid is organized as (nRows x nCols) in BEV coordinates,
/// where rows correspond to Z-axis (leaf travel direction) and
/// columns correspond to X-axis (across leaves).
class FluenceMap {
public:
    FluenceMap() = default;

    /// Extract a 2D fluence map from optimizer weights for a single beam.
    /// @param beam        The beam whose rays define the grid
    /// @param weights     Global optimizer weight vector
    /// @param globalOffset  Global bixel offset for this beam in the weight vector
    static FluenceMap fromBeamWeights(const Beam& beam,
                                     const std::vector<double>& weights,
                                     size_t globalOffset) {
        FluenceMap fm;
        fm.m_bixelWidth = beam.getBixelWidth();

        if (beam.getNumOfRays() == 0) return fm;

        // Determine grid bounds from ray BEV positions
        double xMin =  1e30, xMax = -1e30;
        double zMin =  1e30, zMax = -1e30;
        for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
            const auto& pos = beam.getRay(i)->getRayPosBev();
            xMin = std::min(xMin, pos[0]);
            xMax = std::max(xMax, pos[0]);
            zMin = std::min(zMin, pos[2]);
            zMax = std::max(zMax, pos[2]);
        }

        fm.m_originX = xMin;
        fm.m_originZ = zMin;

        double bw = fm.m_bixelWidth;
        fm.m_nCols = static_cast<int>(std::round((xMax - xMin) / bw)) + 1;
        fm.m_nRows = static_cast<int>(std::round((zMax - zMin) / bw)) + 1;
        fm.m_data.assign(static_cast<size_t>(fm.m_nRows) * fm.m_nCols, 0.0);

        // Fill from weights
        for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
            const auto& pos = beam.getRay(i)->getRayPosBev();
            int col = static_cast<int>(std::round((pos[0] - xMin) / bw));
            int row = static_cast<int>(std::round((pos[2] - zMin) / bw));

            if (row >= 0 && row < fm.m_nRows && col >= 0 && col < fm.m_nCols) {
                size_t weightIdx = globalOffset + i;
                fm.m_data[static_cast<size_t>(row) * fm.m_nCols + col] =
                    (weightIdx < weights.size()) ? weights[weightIdx] : 0.0;
            }
        }

        return fm;
    }

    /// Get fluence value at (row, col).
    double getValue(int row, int col) const {
        if (row < 0 || row >= m_nRows || col < 0 || col >= m_nCols) return 0.0;
        return m_data[static_cast<size_t>(row) * m_nCols + col];
    }

    /// Get a 1D fluence profile for a given row (one per leaf-pair row).
    std::vector<double> getProfile(int row) const {
        std::vector<double> profile(m_nCols, 0.0);
        if (row < 0 || row >= m_nRows) return profile;
        for (int c = 0; c < m_nCols; ++c) {
            profile[c] = m_data[static_cast<size_t>(row) * m_nCols + c];
        }
        return profile;
    }

    /// Map fluence grid rows to MLC leaf pair indices.
    /// For mixed-width MLCs (e.g., Millennium 120: 60x5mm inner + 20x10mm outer per side),
    /// multiple bixel rows may map to one wide leaf pair.
    /// Returns a vector of size nRows, where each element is the leaf pair index.
    std::vector<int> mapToLeafPairs(const MachineGeometry& mlc) const {
        std::vector<int> mapping(m_nRows, -1);
        if (mlc.numLeaves == 0 || mlc.leafWidths.empty()) {
            // No MLC info — 1:1 mapping
            for (int r = 0; r < m_nRows; ++r) mapping[r] = r;
            return mapping;
        }

        // Build leaf pair boundaries.
        // Millennium 120: 60 leaf pairs total (120 leaves = 60 pairs per bank).
        // Organization: 10 outer pairs (10mm) + 40 inner pairs (5mm) + 10 outer pairs (10mm).
        int numPairs = mlc.numLeaves / 2;
        std::vector<double> leafBoundaries; // cumulative Y positions from bottom
        leafBoundaries.push_back(0.0);

        if (mlc.leafWidths.size() == 2) {
            // Mixed width: [inner, outer] = [5mm, 10mm]
            double innerWidth = mlc.leafWidths[0];
            double outerWidth = mlc.leafWidths[1];
            int innerPairs = (mlc.numInnerPairs > 0) ? mlc.numInnerPairs : 40;
            int outerPerSide = (numPairs - innerPairs) / 2;

            // Bottom outer
            for (int i = 0; i < outerPerSide; ++i)
                leafBoundaries.push_back(leafBoundaries.back() + outerWidth);
            // Inner
            for (int i = 0; i < innerPairs; ++i)
                leafBoundaries.push_back(leafBoundaries.back() + innerWidth);
            // Top outer
            for (int i = 0; i < outerPerSide; ++i)
                leafBoundaries.push_back(leafBoundaries.back() + outerWidth);
        } else {
            // Uniform width
            double w = mlc.leafWidths[0];
            for (int i = 0; i < numPairs; ++i)
                leafBoundaries.push_back(leafBoundaries.back() + w);
        }

        double totalLeafSpan = leafBoundaries.back();
        double leafCenterOffset = totalLeafSpan / 2.0;

        // Map each fluence row (by its Z-center position) to the corresponding leaf pair
        for (int r = 0; r < m_nRows; ++r) {
            double zPos = m_originZ + r * m_bixelWidth;
            // Convert to leaf coordinate (0 = bottom leaf edge)
            double leafCoord = zPos + leafCenterOffset;

            for (int lp = 0; lp < numPairs; ++lp) {
                if (leafCoord >= leafBoundaries[lp] && leafCoord < leafBoundaries[lp + 1]) {
                    mapping[r] = lp;
                    break;
                }
            }
        }

        return mapping;
    }

    /// Maximum fluence value across the entire map.
    double getMaxFluence() const {
        if (m_data.empty()) return 0.0;
        return *std::max_element(m_data.begin(), m_data.end());
    }

    int getNumRows() const { return m_nRows; }
    int getNumCols() const { return m_nCols; }
    double getBixelWidth() const { return m_bixelWidth; }
    double getOriginX() const { return m_originX; }
    double getOriginZ() const { return m_originZ; }

    const std::vector<double>& getData() const { return m_data; }

private:
    std::vector<double> m_data;
    int m_nRows = 0;
    int m_nCols = 0;
    double m_bixelWidth = 1.0;
    double m_originX = 0.0; // BEV X position of column 0
    double m_originZ = 0.0; // BEV Z position of row 0
};

} // namespace optirad
