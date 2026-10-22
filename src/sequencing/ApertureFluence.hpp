#pragma once

#include "core/Aperture.hpp"
#include "core/Beam.hpp"
#include <vector>

namespace optirad {

/// Converts an MLC aperture into per-bixel open fractions for a beam's ray grid.
class ApertureFluence {
public:
    /// Fraction of each ray's bixel (bixel width square, BEV x/z) that is open under `aperture`.
    /// bankA/bankB hold the left/right leaf tips per leaf pair; leafPairBoundariesZ has one more entry
    /// than there are leaf pairs. Leaf transmission is ignored. Result is indexed like the beam's rays.
    static std::vector<double> bixelOpenFractions(const Beam& beam,
                                                  const Aperture& aperture,
                                                  const std::vector<double>& leafPairBoundariesZ);

    /// Open fractions of a control point whose leaves keep moving through its angular slice
    /// (integrating continuous motion). Leaf tips follow the piecewise-linear path through the tips of
    /// the previous, current and next control point (positions are those at the slice centres), so the slice
    /// averages the motion from the previous/current midpoint to the current/next midpoint.
    /// `previous`/`next` may be null at the ends of an arc. All apertures must share `leafPairBoundariesZ`.
    static std::vector<double> integratedOpenFractions(const Beam& beam,
                                                       const Aperture* previous,
                                                       const Aperture& current,
                                                       const Aperture* next,
                                                       const std::vector<double>& leafPairBoundariesZ);

    /// Adds d(sum_r rayWeight[r] * integratedOpenFractions[r]) / d(tip) to the tip gradients of the three
    /// apertures. `gradPrevious` / `gradNext` may be null (their share goes to the current aperture,
    /// which is what a missing neighbour means).
    static void accumulateIntegratedTipGradient(const Beam& beam,
                                                const Aperture* previous,
                                                const Aperture& current,
                                                const Aperture* next,
                                                const std::vector<double>& leafPairBoundariesZ,
                                                const std::vector<double>& rayWeight,
                                                Aperture* gradPrevious,
                                                Aperture& gradCurrent,
                                                Aperture* gradNext);
};

} // namespace optirad
