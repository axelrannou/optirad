#pragma once

#include "core/ArcPlan.hpp"
#include "core/Aperture.hpp"
#include <vector>

namespace optirad {

/// A control point with its aperture, MU and the delivery timing derived from machine limits.
struct DeliveredControlPoint {
    ArcControlPoint geometry;
    Aperture aperture;                       // leaf positions in BEV mm; weight = MU
    std::vector<double> leafPairBoundariesZ; // BEV Z edges of the leaf pairs in `aperture` (size = pairs + 1)

    double mu = 0.0;
    double timeSec = 0.0;        // time spent on this slice
    double gantrySpeed = 0.0;    // deg/s actually needed
    double doseRate = 0.0;       // MU/s actually needed
    double leafSpeed = 0.0;      // mm/s, fastest leaf moving into this slice
};

struct ArcSequenceResult {
    std::vector<DeliveredControlPoint> controlPoints;

    double totalMU = 0.0;
    double totalTimeSec = 0.0;
    double maxLeafSpeed = 0.0;       // mm/s
    double meanLeafTravel = 0.0;     // mm per control point, averaged over the leaf tips of the open leaf pairs
    double minGantrySpeed = 0.0;     // deg/s, slowest slice
    int    slowedControlPoints = 0;  // slices where leaf travel or dose rate forced the gantry below its maximum
};

} // namespace optirad
