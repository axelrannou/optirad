#pragma once

#include "core/Machine.hpp"
#include <vector>
#include <cmath>

namespace optirad {

/// Build the full leaf-pair boundary array for a given MLC geometry.
///
/// Returns a vector of size (numPairs + 1) containing the Z-boundaries of each
/// leaf pair in BEV coordinates (mm), centered on isocenter (boundaries[0] < 0,
/// boundaries[numPairs] > 0, midpoint at 0).
///
/// Supported layouts:
///   - Uniform width:    leafWidths.size() == 1  →  all pairs have width leafWidths[0]
///   - Mixed width:      leafWidths.size() == 2  →  numInnerPairs central pairs at
///                       leafWidths[0], remainder split equally on each side at leafWidths[1]
///
/// Examples:
///   Millennium 120: numLeaves=120, leafWidths=[5, 10], numInnerPairs=40
///     → 10 outer (10mm) | 40 inner (5mm) | 10 outer (10mm) = ±200mm
///   HD120:          numLeaves=120, leafWidths=[2.5, 5],  numInnerPairs=32
///     → 14 outer (5mm) | 32 inner (2.5mm) | 14 outer (5mm) = ±110mm
///
/// Returns an empty vector if the MLC geometry is invalid (numLeaves == 0 or
/// leafWidths is empty).
inline std::vector<double> buildMlcLeafBounds(const MachineGeometry& mlc)
{
    int numPairs = (mlc.numLeaves > 0) ? mlc.numLeaves / 2 : 0;
    if (numPairs == 0 || mlc.leafWidths.empty())
        return {};

    std::vector<double> bounds;
    bounds.reserve(static_cast<size_t>(numPairs) + 1);
    bounds.push_back(0.0);

    if (mlc.leafWidths.size() >= 2) {
        // Mixed-width layout: inner (narrow) flanked by outer (wide) pairs
        double innerW = mlc.leafWidths[0];
        double outerW = mlc.leafWidths[1];
        int innerPairs  = (mlc.numInnerPairs > 0) ? mlc.numInnerPairs : 40;
        int outerPerSide = (numPairs - innerPairs) / 2;

        for (int i = 0; i < outerPerSide; ++i)
            bounds.push_back(bounds.back() + outerW);
        for (int i = 0; i < innerPairs; ++i)
            bounds.push_back(bounds.back() + innerW);
        for (int i = 0; i < outerPerSide; ++i)
            bounds.push_back(bounds.back() + outerW);
    } else {
        // Uniform-width layout
        double w = mlc.leafWidths[0];
        for (int i = 0; i < numPairs; ++i)
            bounds.push_back(bounds.back() + w);
    }

    // Center at isocenter (shift so midpoint = 0)
    double halfSpan = bounds.back() / 2.0;
    for (auto& b : bounds)
        b -= halfSpan;

    return bounds;
}

} // namespace optirad
