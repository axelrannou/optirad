#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace optirad {

enum class DeliveryMode { SMLC, VMAT };

/// One gantry arc. Angles are in degrees on the 0..360 scale.
struct ArcSpec {
    double startAngle = 181.0;
    double stopAngle  = 179.0;
    bool   clockwise  = true;   // clockwise = gantry angle increases
    double couchAngle = 0.0;

    /// Angular extent travelled from start to stop (0 < span <= 360).
    double span() const {
        double d = clockwise ? std::fmod(stopAngle - startAngle, 360.0)
                             : std::fmod(startAngle - stopAngle, 360.0);
        if (d < 0.0) d += 360.0;
        return d < 1e-9 ? 360.0 : d;
    }

    /// Gantry angle after travelling `distance` degrees from the start, wrapped to [0, 360).
    double angleAt(double distance) const {
        double a = std::fmod(startAngle + (clockwise ? distance : -distance), 360.0);
        return a < 0.0 ? a + 360.0 : a;
    }
};

struct VmatSettings {
    std::vector<ArcSpec> arcs{ArcSpec{}};
    double optimizationSpacingDeg = 20.0;  // angular width of one fluence-optimization beam
    double controlPointSpacingDeg = 4.0;   // angular width of one delivered control point
};

/// A fluence-optimization beam covering one angular sector of an arc.
struct ArcSector {
    size_t arcIndex = 0;
    double startOffset = 0.0;   // degrees from the arc start
    double width = 0.0;         // degrees
    double centerAngle = 0.0;   // gantry angle at the sector centre
    double couchAngle = 0.0;
};

/// A delivered control point: a short angular slice of an arc with a static aperture.
struct ArcControlPoint {
    size_t arcIndex = 0;
    size_t sectorIndex = 0;     // index into the sector list
    double offset = 0.0;        // degrees from the arc start to the slice centre
    double width = 0.0;         // degrees
    double gantryAngle = 0.0;
    double couchAngle = 0.0;
};

/// Splits arcs into optimization sectors and delivery control points.
/// Sectors are listed arc by arc, in travel order; control points likewise.
class ArcGeometry {
public:
    static std::vector<ArcSector> sectors(const VmatSettings& s) {
        std::vector<ArcSector> out;
        for (size_t a = 0; a < s.arcs.size(); ++a) {
            const ArcSpec& arc = s.arcs[a];
            const double span = arc.span();
            const int n = std::max(1, static_cast<int>(std::lround(span / std::max(s.optimizationSpacingDeg, 1e-6))));
            const double w = span / n;
            for (int i = 0; i < n; ++i) {
                ArcSector sec;
                sec.arcIndex = a;
                sec.startOffset = i * w;
                sec.width = w;
                sec.centerAngle = arc.angleAt((i + 0.5) * w);
                sec.couchAngle = arc.couchAngle;
                out.push_back(sec);
            }
        }
        return out;
    }

    /// Control points, each assigned to the sector that contains its centre.
    static std::vector<ArcControlPoint> controlPoints(const VmatSettings& s) {
        const std::vector<ArcSector> secs = sectors(s);
        std::vector<ArcControlPoint> out;
        size_t firstSector = 0;
        for (size_t a = 0; a < s.arcs.size(); ++a) {
            const ArcSpec& arc = s.arcs[a];
            const double span = arc.span();
            const int n = std::max(1, static_cast<int>(std::lround(span / std::max(s.controlPointSpacingDeg, 1e-6))));
            const double w = span / n;

            size_t sectorsInArc = 0;
            while (firstSector + sectorsInArc < secs.size() && secs[firstSector + sectorsInArc].arcIndex == a)
                ++sectorsInArc;

            for (int j = 0; j < n; ++j) {
                ArcControlPoint cp;
                cp.arcIndex = a;
                cp.offset = (j + 0.5) * w;
                cp.width = w;
                cp.gantryAngle = arc.angleAt(cp.offset);
                cp.couchAngle = arc.couchAngle;
                size_t local = 0;
                while (local + 1 < sectorsInArc &&
                       cp.offset >= secs[firstSector + local].startOffset + secs[firstSector + local].width)
                    ++local;
                cp.sectorIndex = firstSector + local;
                out.push_back(cp);
            }
            firstSector += sectorsInArc;
        }
        return out;
    }
};

} // namespace optirad
