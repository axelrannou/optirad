#pragma once

namespace optirad {

/**
 * Options for dose influence matrix (Dij) computation.
 *
 * Controls memory usage and parallelism during dose calculation.
 * Inspired by matRad's batched dose calculation with thresholding.
 */
struct DoseCalcOptions {
    /// Discard dose entries below this absolute value (Gy).
    /// Removes noise from far-field kernel tails. Default 1e-6.
    double absoluteThreshold = 1e-6;

    /// Discard dose entries below this fraction of the bixel's max dose.
    /// Applied per-bixel after computing all voxel doses for that bixel.
    /// 1e-4 = keep only values >= 0.01% of max. Scatter tails of many bixels add up, so
    /// 1% removes ~15% of the dose at depth (validated on the water phantom).
    double relativeThreshold = 1e-4;

    /// Skip voxels that belong only to the external/body structure (saves memory).
    /// Dose there is only available on a coarse lattice; use calculateDoseDirect() for the full cube.
    bool excludeExternal = true;

    /// With excludeExternal, keep every n-th voxel (per axis) of the body so dose-wide terms such as the
    /// hotspot penalty still see the dose outside the structures. 1 = keep none. Changing it alters the Dij
    /// content, so bump DijSerializer::kEngineVersion.
    int externalSampleStep = 3;

    /// Number of OpenMP threads. 0 = use all available (default).
    int numThreads = 0;
};

} // namespace optirad
