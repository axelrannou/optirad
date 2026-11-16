#pragma once

#include "core/ArcPlan.hpp"
#include "core/ArcSequence.hpp"
#include "core/Machine.hpp"
#include "core/Stf.hpp"
#include <vector>

namespace optirad {

struct ArcSequencerOptions {
    /// A leaf pair is open where the sector fluence exceeds this fraction of the sector's peak.
    double conformalThreshold = 0.1;
    double leafPositionResolution = 0.5; // mm
};

/// Control point apertures for VMAT direct aperture optimization.
///
/// The conformal seed opens every leaf pair over the target projection of each sector's beam (the rays of
/// `sectorStf`), interpolated between sector centres so the leaves glide along the arc. The modulation and the MU are
/// left to the DAO of ArcPipeline.
class ArcSequencer {
public:
    /// `sectorStf` beams must be the sector centres of `settings` in the order produced by ArcGeometry::sectors().
    static ArcSequenceResult conformalSeed(const Stf& sectorStf,
                                           const VmatSettings& settings,
                                           const Machine& machine,
                                           const ArcSequencerOptions& options = {});

    /// Resamples a sequence onto other control points of the same arcs (progressive resolution): leaf tips are
    /// interpolated linearly in the arc angle between the surrounding control points, MU per degree likewise
    /// (constant beyond the first/last control point of an arc). Timing is recomputed with `limits`.
    static ArcSequenceResult resample(const ArcSequenceResult& source,
                                      const std::vector<ArcControlPoint>& geometry,
                                      const MachineConstraints& limits);

    /// Recompute timing and totals from the control points' MU and aperture positions.
    /// Gantry speed, dose rate and leaf speed are bounded by the machine constraints; a slice is slowed
    /// (the gantry waits) when leaf travel or the MU rate cannot keep up with the maximum gantry speed.
    static void computeTiming(ArcSequenceResult& result, const MachineConstraints& limits);
};

} // namespace optirad
