#pragma once

#include "core/ArcPlan.hpp"
#include "core/ArcSequence.hpp"
#include "core/Machine.hpp"
#include "core/Stf.hpp"
#include <vector>

namespace optirad {

struct ArcSequencerOptions {
    enum class Method {
        Auto,      // sliding window when the MLC can interdigitate and every sector has >= 2 control points
        Sweep,     // sliding window (throws if not applicable)
        LevelSet,  // intensity-level decomposition, one segment per control point
        Conformal  // apertures following the target outline in the beam's eye view; leaf optimization adds the modulation
    };
    Method method = Method::Auto;

    /// Sliding window: the sweep speed and the MU per control point are chosen together to minimize the MU of a
    /// sector: the leaves travel at the machine's leaf speed over the (possibly slowed) slice time, and the MU
    /// that remain are delivered while they move. This caps the share of a sector's time spent sweeping.
    double maxSweepFraction = 0.6;

    /// Sliding window: a leaf pair delivering a very modulated row needs MU in proportion to the sum of its upward
    /// fluence steps, and that sets the MU of the whole sector. Rows whose steps exceed this multiple of the
    /// sector's peak fluence are smoothed (mass preserving) until they fit. 0 = no smoothing.
    double maxRiseRatio = 1.25;

    /// Sweep: leaf tips at the slice centres are refined (projected gradient, leaf speed respected) so the
    /// fluence delivered by the moving leaves matches the optimized fluence of each sector. 0 = use the
    /// sampled sweep trajectories as they are.
    int fitIterations = 80;

    /// Conformal: a leaf pair is open where the sector fluence exceeds this fraction of the sector's peak.
    double conformalThreshold = 0.1;

    int numLevels = 15;                  // LevelSet: upper bound; limited per sector by its number of control points
    double minSegmentMU = 0.0;
    double leafPositionResolution = 0.5; // mm
    bool enforceLeafSpeed = true;        // bend leaf trajectories so they are reachable at full gantry speed
};

/// Arc sequencing: turns the fluence of the optimization beams (one per arc sector) into one aperture per
/// control point.
///
/// Sweep (default): every leaf pair delivers its sector's fluence profile with a sliding window. Consecutive
/// sectors sweep in opposite directions and meet at a common leaf position, so the leaves move continuously
/// along the whole arc. Control point apertures hold the leaf tips at the slice centres; the leaves are
/// meant to keep moving between them (ApertureFluence::integratedOpenFractions).
///
/// Conformal: every control point opens each leaf pair over the extent of the sector's fluence (linearly
/// interpolated between sector centres), so the leaves stay near the target outline and move little. The
/// modulation is left to the leaf and MU optimization of ArcPipeline; MU are scaled to the sector fluence.
///
/// LevelSet: each sector's fluence is decomposed into at most as many segments as it has control
/// points; segments are spread over the sector's control points and traversed in alternating order
/// (large to small, then small to large) so the leaves move continuously across sector boundaries.
class ArcSequencer {
public:
    /// `weights` are the optimized bixel weights of `fmoStf`, whose beams must be the sector centres of
    /// `settings` in the order produced by ArcGeometry::sectors().
    static ArcSequenceResult sequence(const std::vector<double>& weights,
                                      const Stf& fmoStf,
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
