#pragma once

#include <vector>
#include <string>

namespace optirad {

/// A single MLC aperture segment for step-and-shoot IMRT delivery.
/// Each segment defines leaf bank positions (at isocenter plane) and a MU weight.
struct Aperture {
    /// Left leaf bank positions (mm at isocenter plane, one per leaf pair).
    /// Positive values = leaf extends from left towards center.
    std::vector<double> bankA;

    /// Right leaf bank positions (mm at isocenter plane, one per leaf pair).
    /// Positive values = leaf extends from right towards center.
    std::vector<double> bankB;

    /// Monitor Units for this segment.
    double weight = 0.0;

    /// Segment index within the beam.
    int segmentIndex = 0;
};

/// Result of leaf sequencing for a single beam.
struct LeafSequenceResult {
    /// Ordered list of aperture segments.
    std::vector<Aperture> segments;

    /// Beam index in the Stf.
    size_t beamIndex = 0;

    /// Total MU across all segments.
    double totalMU = 0.0;

    /// Fluence fidelity: correlation between deliverable and optimal fluence.
    double fluenceFidelity = 0.0;

    /// Quantized fluence: D_0 matrix stored in row-major order (nRows x nCols).
    /// Used for exact deliverable weight reconstruction: w = D_0(idx) / numLevels * calFac.
    std::vector<int> quantizedFluence;

    /// Calibration factor (max fluence before quantization).
    double calFac = 0.0;

    /// Number of quantization levels used.
    int numLevels = 0;

    /// Z boundaries (BEV mm) for each leaf pair in the field.
    /// Size = numLeafPairsInField + 1.  leafPairBoundariesZ[i] is the bottom
    /// edge of leaf pair i, leafPairBoundariesZ[i+1] is its top edge.
    /// Empty when MLC info is unavailable (callers fall back to bixel width).
    std::vector<double> leafPairBoundariesZ;

    /// Deliverable fluence at leaf-pair resolution (row-major, numFieldLP × leafPairFluenceCols).
    /// Used by BevView to render a texture aligned with physical leaf pair widths.
    std::vector<double> leafPairFluence;

    /// Number of columns in leafPairFluence.
    int leafPairFluenceCols = 0;

    /// BEV X origin (mm) for the first column of the fluence grid.
    double originX = 0.0;
};

/// Configuration options for the leaf sequencer.
struct LeafSequencerOptions {
    /// Number of discrete intensity levels for fluence quantization.
    int numLevels = 15;

    /// Minimum MU a segment must have to survive.
    /// Segments below this threshold are merged with their neighbor.
    /// 0 = disabled (keep all segments).
    double minSegmentMU = 0.0;

    /// Snap leaf positions to this resolution (mm).
    double leafPositionResolution = 0.5;
};

} // namespace optirad
