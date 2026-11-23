#include "ImportedPlanConverter.hpp"
#include "core/Machine.hpp"
#include "core/Plan.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "sequencing/ArcSequencer.hpp"
#include "sequencing/MlcLeafBounds.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace optirad {

namespace {

/// Gantry travel (degrees, [0, 360)) from `from` to `to` in the given rotation direction.
double travel(double from, double to, bool clockwise) {
    double d = std::fmod(clockwise ? to - from : from - to, 360.0);
    if (d < 0.0) d += 360.0;
    return d > 360.0 - 1e-6 ? 0.0 : d;
}

double arcTravel(const ImportedBeam& beam) {
    double sum = 0.0;
    for (size_t k = 0; k + 1 < beam.controlPoints.size(); ++k)
        sum += travel(beam.controlPoints[k].gantryAngle, beam.controlPoints[k + 1].gantryAngle, beam.clockwise);
    return sum;
}

ImportedControlPoint midpoint(const ImportedControlPoint& a, const ImportedControlPoint& b) {
    ImportedControlPoint m = a;
    if (a.bankA.size() == b.bankA.size() && a.bankB.size() == b.bankB.size()) {
        for (size_t i = 0; i < m.bankA.size(); ++i) m.bankA[i] = 0.5 * (a.bankA[i] + b.bankA[i]);
        for (size_t i = 0; i < m.bankB.size(); ++i) m.bankB[i] = 0.5 * (a.bankB[i] + b.bankB[i]);
    }
    for (int i = 0; i < 2; ++i) {
        m.jawX[i] = 0.5 * (a.jawX[i] + b.jawX[i]);
        m.jawY[i] = 0.5 * (a.jawY[i] + b.jawY[i]);
    }
    return m;
}

} // namespace

Aperture ImportedPlanConverter::collimatedAperture(const ImportedControlPoint& cp,
                                                   const std::vector<double>& leafBoundaries) {
    const size_t pairs = leafBoundaries.size() >= 2 ? leafBoundaries.size() - 1 : 0;
    const bool hasLeaves = cp.bankA.size() == pairs && cp.bankB.size() == pairs;
    Aperture ap;
    ap.bankA.resize(pairs);
    ap.bankB.resize(pairs);
    for (size_t i = 0; i < pairs; ++i) {
        const double a0 = hasLeaves ? cp.bankA[i] : cp.jawX[0];
        const double b0 = hasLeaves ? cp.bankB[i] : cp.jawX[1];
        double a = std::max(a0, cp.jawX[0]);
        double b = std::min(b0, cp.jawX[1]);
        const bool shielded = leafBoundaries[i + 1] <= cp.jawY[0] || leafBoundaries[i] >= cp.jawY[1];
        if (shielded || b <= a) a = b = std::clamp(0.5 * (a0 + b0), cp.jawX[0], cp.jawX[1]);
        ap.bankA[i] = a;
        ap.bankB[i] = b;
    }
    return ap;
}

VmatSettings ImportedPlanConverter::toVmatSettings(const ImportedPlan& imported) {
    VmatSettings s;
    s.arcs.clear();
    double cpSpacing = std::numeric_limits<double>::infinity();
    for (const auto& beam : imported.beams) {
        if (!beam.isArc) continue;
        ArcSpec arc;
        arc.startAngle = beam.controlPoints.front().gantryAngle;
        arc.stopAngle = beam.controlPoints.back().gantryAngle;
        arc.clockwise = beam.clockwise;
        arc.couchAngle = beam.couchAngle;
        arc.collimatorAngle = beam.collimatorAngle;
        s.arcs.push_back(arc);
        cpSpacing = std::min(cpSpacing, arcTravel(beam) / static_cast<double>(beam.controlPoints.size() - 1));
    }
    if (s.arcs.empty() || !(cpSpacing > 1e-6)) return s;
    s.controlPointSpacingDeg = cpSpacing;
    s.optimizationSpacingDeg = cpSpacing * static_cast<double>(std::max(1L, std::lround(10.0 / cpSpacing)));
    return s;
}

void ImportedPlanConverter::applyToPlanConfig(const ImportedPlan& imported, PlanConfig& config) {
    if (imported.numFractions > 0) config.numFractions = imported.numFractions;
    config.isocenter = imported.isocenter();
    if (imported.prescribedDose > 0.0) config.prescribedDose = imported.prescribedDose;

    if (imported.hasArcs()) {
        config.deliveryMode = DeliveryMode::VMAT;
        config.vmat = toVmatSettings(imported);
        return;
    }
    config.deliveryMode = DeliveryMode::SMLC;
    config.gantryAngles.clear();
    config.couchAngles.clear();
    config.collimatorAngles.clear();
    for (const auto& beam : imported.beams) {
        config.gantryAngles.push_back(beam.controlPoints.front().gantryAngle);
        config.couchAngles.push_back(beam.couchAngle);
        config.collimatorAngles.push_back(beam.collimatorAngle);
    }
}

ArcSequenceResult ImportedPlanConverter::toArcSequence(const ImportedPlan& imported, const Plan& plan) {
    VmatSettings settings = toVmatSettings(imported);
    if (plan.isVmat() && plan.getVmatSettings().arcs.size() == settings.arcs.size() &&
        plan.getVmatSettings().optimizationSpacingDeg > 0.0)
        settings.optimizationSpacingDeg = plan.getVmatSettings().optimizationSpacingDeg;
    const std::vector<ArcSector> sectors = ArcGeometry::sectors(settings);
    const double muPerWeight = ArcPipeline::machineMuPerWeight(plan);

    ArcSequenceResult result;
    size_t arcIndex = 0, firstSector = 0;
    for (const auto& beam : imported.beams) {
        if (!beam.isArc) continue;
        const ArcSpec& arc = settings.arcs[arcIndex];
        size_t sectorsInArc = 0;
        while (firstSector + sectorsInArc < sectors.size() && sectors[firstSector + sectorsInArc].arcIndex == arcIndex)
            ++sectorsInArc;
        const double sectorWidth = sectorsInArc > 0 ? sectors[firstSector].width : arc.span();

        const auto& cps = beam.controlPoints;
        double offset = 0.0;
        double pendingMu = 0.0;   // MU delivered while the gantry stands still, given to the next slice
        for (size_t j = 0; j + 1 < cps.size(); ++j) {
            const double width = travel(cps[j].gantryAngle, cps[j + 1].gantryAngle, beam.clockwise);
            const double mu = std::max(0.0, cps[j + 1].cumulativeWeight - cps[j].cumulativeWeight) *
                              beam.metersetMU / muPerWeight;
            if (width < 1e-6) { pendingMu += mu; continue; }

            DeliveredControlPoint dcp;
            dcp.geometry.arcIndex = arcIndex;
            dcp.geometry.offset = offset + 0.5 * width;
            dcp.geometry.width = width;
            dcp.geometry.gantryAngle = arc.angleAt(dcp.geometry.offset);
            dcp.geometry.couchAngle = arc.couchAngle;
            dcp.geometry.collimatorAngle = arc.collimatorAngle;
            const size_t local = std::min(sectorsInArc > 0 ? sectorsInArc - 1 : 0,
                                          static_cast<size_t>(dcp.geometry.offset / sectorWidth));
            dcp.geometry.sectorIndex = firstSector + local;

            dcp.aperture = collimatedAperture(midpoint(cps[j], cps[j + 1]), beam.leafBoundaries);
            dcp.mu = mu + pendingMu;
            pendingMu = 0.0;
            dcp.aperture.weight = dcp.mu;
            dcp.aperture.segmentIndex = static_cast<int>(result.controlPoints.size());
            dcp.leafPairBoundariesZ = beam.leafBoundaries;
            result.controlPoints.push_back(std::move(dcp));
            offset += width;
        }
        if (pendingMu > 0.0 && !result.controlPoints.empty()) {
            result.controlPoints.back().mu += pendingMu;
            result.controlPoints.back().aperture.weight = result.controlPoints.back().mu;
        }
        firstSector += sectorsInArc;
        ++arcIndex;
    }
    ArcSequencer::computeTiming(result, ArcPipeline::deliveryConstraints(plan));
    return result;
}

std::vector<std::string> ImportedPlanConverter::arcPlanMismatches(const ImportedPlan& imported, const Plan& plan) {
    std::vector<std::string> m;
    if (!imported.hasArcs()) {
        m.push_back("the RTPLAN has no arc beams");
        return m;
    }
    if (!plan.isVmat()) {
        m.push_back("the plan is not a VMAT plan");
        return m;
    }
    auto angleDiff = [](double a, double b) {
        const double d = std::fmod(std::abs(a - b), 360.0);
        return std::min(d, 360.0 - d);
    };
    constexpr double tolDeg = 1e-3;
    const auto& want = toVmatSettings(imported).arcs;
    const auto& have = plan.getVmatSettings().arcs;
    if (want.size() != have.size()) {
        m.push_back("the plan has " + std::to_string(have.size()) + " arc(s), the RTPLAN " +
                    std::to_string(want.size()));
    } else {
        for (size_t a = 0; a < want.size(); ++a) {
            const auto& w = want[a];
            const auto& h = have[a];
            if (angleDiff(w.startAngle, h.startAngle) > tolDeg || angleDiff(w.stopAngle, h.stopAngle) > tolDeg ||
                w.clockwise != h.clockwise || std::abs(w.span() - h.span()) > tolDeg)
                m.push_back("arc " + std::to_string(a + 1) + " gantry range or direction differs from the RTPLAN");
            if (angleDiff(w.couchAngle, h.couchAngle) > tolDeg)
                m.push_back("arc " + std::to_string(a + 1) + " couch angle differs from the RTPLAN");
            if (angleDiff(w.collimatorAngle, h.collimatorAngle) > tolDeg)
                m.push_back("arc " + std::to_string(a + 1) + " collimator angle differs from the RTPLAN");
        }
    }
    if (imported.numFractions > 0 && plan.getNumOfFractions() != imported.numFractions)
        m.push_back("the plan has " + std::to_string(plan.getNumOfFractions()) + " fraction(s), the RTPLAN " +
                    std::to_string(imported.numFractions));
    const auto& isoCenters = plan.getStfProperties().isoCenters;
    if (const auto iso = imported.isocenter()) {
        if (isoCenters.empty()) {
            m.push_back("the plan has no isocenter");
        } else {
            const auto& p = isoCenters.front();
            const double d = std::sqrt(std::pow(p[0] - (*iso)[0], 2) + std::pow(p[1] - (*iso)[1], 2) +
                                       std::pow(p[2] - (*iso)[2], 2));
            if (d > 0.1) m.push_back("the plan isocenter is " + std::to_string(d) + " mm from the RTPLAN's");
        }
    }
    return m;
}

std::vector<LeafSequenceResult> ImportedPlanConverter::toLeafSequences(const ImportedPlan& imported,
                                                                      const Plan& plan) {
    const double muPerWeight = ArcPipeline::machineMuPerWeight(plan);
    std::vector<LeafSequenceResult> out;
    for (const auto& beam : imported.beams) {
        if (beam.isArc) continue;
        LeafSequenceResult seq;
        seq.beamIndex = out.size();
        seq.leafPairBoundariesZ = beam.leafBoundaries;
        const auto& cps = beam.controlPoints;
        for (size_t k = 0; k + 1 < cps.size(); ++k) {
            const double dw = cps[k + 1].cumulativeWeight - cps[k].cumulativeWeight;
            if (dw <= 1e-9) continue;
            Aperture ap = collimatedAperture(cps[k], beam.leafBoundaries);
            ap.weight = dw * beam.metersetMU / muPerWeight;
            ap.segmentIndex = static_cast<int>(seq.segments.size());
            seq.totalMU += ap.weight;
            seq.segments.push_back(std::move(ap));
        }
        out.push_back(std::move(seq));
    }
    return out;
}

std::vector<std::string> ImportedPlanConverter::validate(const ImportedPlan& imported, const Machine& machine) {
    std::vector<std::string> w;
    if (imported.patientPosition != "HFS")
        w.push_back("patient position " + imported.patientPosition + ": leaf directions assume HFS");
    if (imported.hasMultipleIsocenters())
        w.push_back("several isocenters: optirad uses the first beam's isocenter for every beam");
    bool hasStatic = false;
    for (const auto& b : imported.beams) hasStatic |= !b.isArc;
    if (imported.hasArcs() && hasStatic)
        w.push_back("plan mixes arcs and static beams: the static beams are ignored in VMAT mode");

    const auto& geo = machine.getGeometry();
    const std::vector<double> machineBounds = buildMlcLeafBounds(geo);
    for (const auto& b : imported.beams) {
        const std::string tag = "beam '" + b.name + "': ";
        if (!b.radiationType.empty() && b.radiationType != "PHOTON")
            w.push_back(tag + b.radiationType + " beams are not supported");
        if (b.fff)
            w.push_back(tag + "flattening-filter-free beam, the machine model is flattened");
        double normalizedCollimator = std::fmod(b.collimatorAngle, 360.0);
        if (normalizedCollimator > 180.0) normalizedCollimator -= 360.0;
        if (normalizedCollimator < -180.0) normalizedCollimator += 360.0;
        if (normalizedCollimator < geo.collimatorRange[0] - 1e-3 ||
            normalizedCollimator > geo.collimatorRange[1] + 1e-3)
            w.push_back(tag + "collimator " + std::to_string(b.collimatorAngle) +
                        " deg is outside the machine's supported range");
        if (std::abs(b.couchAngle) > 1e-3)
            w.push_back(tag + "couch " + std::to_string(b.couchAngle) + " deg: IEC couch convention not verified");
        if (machine.getSAD() > 0.0 && std::abs(b.sad - machine.getSAD()) > 0.5)
            w.push_back(tag + "SAD " + std::to_string(b.sad) + " mm differs from the machine's " +
                        std::to_string(machine.getSAD()) + " mm");
        if (geo.doseRateMUPerMin > 0.0 && b.doseRateMUPerMin > geo.doseRateMUPerMin + 1e-6)
            w.push_back(tag + "dose rate " + std::to_string(b.doseRateMUPerMin) + " MU/min exceeds the machine's " +
                        std::to_string(geo.doseRateMUPerMin) + " MU/min");
        if (!b.machineName.empty() && b.machineName != machine.getName())
            w.push_back(tag + "planned on '" + b.machineName + "', recomputed with '" + machine.getName() + "'");

        if (b.leafBoundaries.empty() || b.controlPoints.front().bankA.empty()) {
            w.push_back(tag + "no MLCX leaves, the field is shaped by the jaws only");
        } else if (b.leafBoundaries.size() != machineBounds.size()) {
            w.push_back(tag + std::to_string(b.leafBoundaries.size() - 1) + " leaf pairs, the machine has " +
                        std::to_string(machineBounds.empty() ? 0 : machineBounds.size() - 1));
        } else {
            double maxDiff = 0.0;
            for (size_t i = 0; i < machineBounds.size(); ++i)
                maxDiff = std::max(maxDiff, std::abs(b.leafBoundaries[i] - machineBounds[i]));
            if (maxDiff > 0.1)
                w.push_back(tag + "leaf boundaries differ from the machine's by up to " + std::to_string(maxDiff) + " mm");
        }

        int partial = 0;
        for (const auto& cp : b.controlPoints)
            for (size_t i = 0; i + 1 < b.leafBoundaries.size(); ++i) {
                const double lo = b.leafBoundaries[i], hi = b.leafBoundaries[i + 1];
                if ((cp.jawY[0] > lo && cp.jawY[0] < hi) || (cp.jawY[1] > lo && cp.jawY[1] < hi)) { ++partial; break; }
            }
        if (partial > 0)
            w.push_back(tag + "a Y jaw cuts through a leaf pair at " + std::to_string(partial) +
                        " control point(s); the partly covered pair is kept fully open");

        if (!b.isArc) {
            bool dynamic = false;
            for (size_t k = 0; k + 1 < b.controlPoints.size(); ++k)
                if (b.controlPoints[k + 1].cumulativeWeight > b.controlPoints[k].cumulativeWeight &&
                    (b.controlPoints[k].bankA != b.controlPoints[k + 1].bankA ||
                     b.controlPoints[k].bankB != b.controlPoints[k + 1].bankB))
                    dynamic = true;
            if (dynamic)
                w.push_back(tag + "leaves move while the beam is on (sliding window); each segment uses its "
                                  "starting aperture");
        }
    }
    return w;
}

} // namespace optirad
