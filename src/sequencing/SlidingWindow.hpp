#pragma once

#include <vector>

namespace optirad {

/// One leaf pair delivering a 1D fluence profile with a unidirectional sliding window.
/// The leading leaf waits where the profile falls and the trailing leaf waits where it rises, so every
/// pixel stays exposed for (profile / MU rate) while both leaves travel at most `speed` per time slot.
/// Coordinates are along the sweep direction (increasing); callers mirror the profile for the opposite sweep.
struct SlidingWindowRow {
    std::vector<double> profile;   // MU delivered to each pixel
    double firstCenter = 0.0;      // position of pixel 0's centre
    double pixel = 1.0;            // pixel width
    double start = 0.0;            // both leaves start here (<= left edge of the first exposed pixel)
    double end = 0.0;              // both leaves finish here (>= right edge of the last exposed pixel)
};

class SlidingWindow {
public:
    /// Total MU a sweep must deliver while the window is moving (sum of the upward steps of the profile).
    /// The leaf pair needs at least this much MU in addition to the time it spends travelling.
    static double totalRise(const std::vector<double>& profile);

    /// MU per slot needed so the row fits in `slots` time slots with leaf speed `speed` (distance per slot).
    /// 0 for an empty row; negative if the sweep alone does not fit in `slots`.
    static double requiredMuPerSlot(const SlidingWindowRow& row, double speed, int slots);

    /// Leaf tips at the centre of each slot (slot k covers time [k, k+1)) when `muPerSlot` MU are delivered
    /// per slot. `trailing` is the leaf that starts the exposure last, `leading` the one that starts first.
    static void sampleTips(const SlidingWindowRow& row, double speed, double muPerSlot, int slots,
                           std::vector<double>& trailing, std::vector<double>& leading);
};

} // namespace optirad
