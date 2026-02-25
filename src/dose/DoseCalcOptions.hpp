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
    /// 0.01 = keep only values >= 1% of max. Default 0.01.
    double relativeThreshold = 0.01;

    /// Number of OpenMP threads. 0 = use all available (default).
    int numThreads = 0;
};

} // namespace optirad
