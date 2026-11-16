#include "ArcSequencer.hpp"

#include "core/FluenceMap.hpp"
#include "ApertureFluence.hpp"
#include "MlcLeafBounds.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace optirad {

namespace {

// Largest leaf-tip displacement between two apertures, matching leaf pairs by their Z position.
// Closed pairs still have physical tip positions and contribute to travel.
double maxLeafTravel(const DeliveredControlPoint& a, const DeliveredControlPoint& b) {
    const auto& za = a.leafPairBoundariesZ;
    const auto& zb = b.leafPairBoundariesZ;
    const size_t na = std::min(a.aperture.bankA.size(), za.empty() ? size_t(0) : za.size() - 1);
    const size_t nb = std::min(b.aperture.bankA.size(), zb.empty() ? size_t(0) : zb.size() - 1);

    double travel = 0.0;
    for (size_t i = 0; i < nb; ++i) {
        const double bA = b.aperture.bankA[i], bB = b.aperture.bankB[i];
        const double zc = 0.5 * (zb[i] + zb[i + 1]);
        for (size_t j = 0; j < na; ++j) {
            if (zc < za[j] || zc >= za[j + 1]) continue;
            const double aA = a.aperture.bankA[j], aB = a.aperture.bankB[j];
            travel = std::max(travel, std::max(std::abs(bA - aA), std::abs(bB - aB)));
            break;
        }
    }
    return travel;
}

// One sector's fluence resampled to the MLC leaf pairs (rows) and its per-row exposed extent.
struct SectorRows {
    std::vector<std::vector<double>> profile;   // [pair][column]
    std::vector<double> left, right;            // exposed extent per pair (x edges); 0,0 when the pair is empty
    double originX = 0.0, pixel = 1.0;
    int columns = 0;
};

SectorRows resampleToLeafPairs(const FluenceMap& fluence, const std::vector<double>& leafBounds) {
    SectorRows rows;
    const size_t pairs = leafBounds.size() - 1;
    rows.columns = fluence.getNumCols();
    rows.originX = fluence.getOriginX();
    rows.pixel = fluence.getBixelWidth();
    rows.profile.assign(pairs, std::vector<double>(static_cast<size_t>(std::max(rows.columns, 0)), 0.0));
    rows.left.assign(pairs, 0.0);
    rows.right.assign(pairs, 0.0);
    if (rows.columns <= 0 || fluence.getNumRows() <= 0) return rows;

    const double bw = rows.pixel;
    for (size_t pair = 0; pair < pairs; ++pair) {
        const double width = leafBounds[pair + 1] - leafBounds[pair];
        if (width <= 0.0) continue;
        for (int r = 0; r < fluence.getNumRows(); ++r) {
            const double rowBottom = fluence.getOriginZ() + r * bw - 0.5 * bw;
            const double overlap = std::min(leafBounds[pair + 1], rowBottom + bw) - std::max(leafBounds[pair], rowBottom);
            if (overlap <= 0.0) continue;
            const double share = overlap / width;
            for (int c = 0; c < rows.columns; ++c) rows.profile[pair][static_cast<size_t>(c)] += share * fluence.getValue(r, c);
        }
    }

    // Ignore negligible fluence so it does not widen the aperture.
    double peak = 0.0;
    for (const auto& row : rows.profile) for (double v : row) peak = std::max(peak, v);
    const double floor = 0.005 * peak;
    for (size_t pair = 0; pair < pairs; ++pair) {
        int first = -1, last = -1;
        for (int c = 0; c < rows.columns; ++c) {
            double& v = rows.profile[pair][static_cast<size_t>(c)];
            if (v <= floor) { v = 0.0; continue; }
            if (first < 0) first = c;
            last = c;
        }
        if (first < 0) continue;
        rows.left[pair] = rows.originX + (first - 0.5) * bw;
        rows.right[pair] = rows.originX + (last + 0.5) * bw;
    }
    return rows;
}

// Conformal seed: each leaf pair is open over the extent of the sector fluence above a threshold; the extents are
// interpolated between sector centres so the leaves glide along the arc. MU per control point is the least-squares
// scale that matches the sector's fluence.
void conformalSequence(ArcSequenceResult& result,
                       const std::vector<ArcSector>& sectors,
                       const std::vector<std::vector<size_t>>& cpsOfSector,
                       const std::vector<SectorRows>& rows,
                       const std::vector<double>& leafBounds,
                       const Stf& fmoStf,
                       const std::vector<double>& weights,
                       const std::vector<size_t>& offset,
                       const ArcSequencerOptions& options) {
    const size_t pairs = leafBounds.size() - 1;
    struct Extent { double left = 0.0, right = 0.0; };
    std::vector<std::vector<Extent>> extent(sectors.size(), std::vector<Extent>(pairs));
    for (size_t s = 0; s < sectors.size(); ++s) {
        const SectorRows& sector = rows[s];
        double peak = 0.0;
        for (const auto& row : sector.profile) for (double v : row) peak = std::max(peak, v);
        const double threshold = options.conformalThreshold * peak;
        std::vector<char> open(pairs, 0);
        for (size_t p = 0; p < pairs; ++p) {
            int first = -1, last = -1;
            for (int c = 0; c < sector.columns; ++c) {
                if (sector.profile[p][static_cast<size_t>(c)] <= threshold) continue;
                if (first < 0) first = c;
                last = c;
            }
            if (first < 0) continue;
            open[p] = 1;
            extent[s][p] = {sector.originX + (first - 0.5) * sector.pixel, sector.originX + (last + 0.5) * sector.pixel};
        }
        // Closed pairs tuck in next to the closest open pair so they do not have to travel far.
        for (size_t p = 0; p < pairs; ++p) {
            if (open[p]) continue;
            for (size_t d = 1; d < pairs; ++d) {
                const long below = static_cast<long>(p) - static_cast<long>(d);
                const size_t above = p + d;
                const long near = (below >= 0 && open[static_cast<size_t>(below)]) ? below :
                    (above < pairs && open[above] ? static_cast<long>(above) : -1);
                if (near < 0) continue;
                const double middle = 0.5 * (extent[s][static_cast<size_t>(near)].left + extent[s][static_cast<size_t>(near)].right);
                extent[s][p] = {middle, middle};
                break;
            }
        }
    }

    const double resolution = options.leafPositionResolution;
    auto snap = [resolution](double v) { return resolution > 0.0 ? std::round(v / resolution) * resolution : v; };
    size_t arcStart = 0;
    while (arcStart < sectors.size()) {
        size_t arcEnd = arcStart + 1;
        while (arcEnd < sectors.size() && sectors[arcEnd].arcIndex == sectors[arcStart].arcIndex) ++arcEnd;
        for (size_t s = arcStart; s < arcEnd; ++s) {
            for (size_t j : cpsOfSector[s]) {
                DeliveredControlPoint& cp = result.controlPoints[j];
                const double position = cp.geometry.offset;
                size_t lower = s;
                auto centre = [&](size_t k) { return sectors[k].startOffset + 0.5 * sectors[k].width; };
                if (position < centre(s) && s > arcStart) lower = s - 1;
                const size_t upper = std::min(lower + 1, arcEnd - 1);
                const double span = centre(upper) - centre(lower);
                const double t = span > 0.0 ? std::clamp((position - centre(lower)) / span, 0.0, 1.0) : 0.0;
                cp.leafPairBoundariesZ = leafBounds;
                cp.aperture.bankA.assign(pairs, 0.0);
                cp.aperture.bankB.assign(pairs, 0.0);
                for (size_t p = 0; p < pairs; ++p) {
                    cp.aperture.bankA[p] = snap((1.0 - t) * extent[lower][p].left + t * extent[upper][p].left);
                    cp.aperture.bankB[p] = snap((1.0 - t) * extent[lower][p].right + t * extent[upper][p].right);
                    if (cp.aperture.bankA[p] > cp.aperture.bankB[p]) cp.aperture.bankA[p] = cp.aperture.bankB[p];
                }
                cp.aperture.segmentIndex = static_cast<int>(j);
            }
        }
        arcStart = arcEnd;
    }

    // Per-sector MU: least-squares scale of the delivered unit fluence onto the optimized fluence.
    for (size_t s = 0; s < sectors.size(); ++s) {
        const auto& slots = cpsOfSector[s];
        const Beam& beam = *fmoStf.getBeam(s);
        std::vector<double> unit(beam.getNumOfRays(), 0.0);
        for (size_t j : slots) {
            const Aperture* before = (j > 0 && result.controlPoints[j - 1].geometry.arcIndex == result.controlPoints[j].geometry.arcIndex) ?
                &result.controlPoints[j - 1].aperture : nullptr;
            const Aperture* after = (j + 1 < result.controlPoints.size() && result.controlPoints[j + 1].geometry.arcIndex == result.controlPoints[j].geometry.arcIndex) ?
                &result.controlPoints[j + 1].aperture : nullptr;
            const auto open = ApertureFluence::integratedOpenFractions(beam, before, result.controlPoints[j].aperture, after,
                                                                       result.controlPoints[j].leafPairBoundariesZ);
            for (size_t r = 0; r < open.size(); ++r) unit[r] += open[r];
        }
        double cross = 0.0, norm = 0.0;
        for (size_t r = 0; r < unit.size(); ++r) {
            cross += unit[r] * weights[offset[s] + r];
            norm += unit[r] * unit[r];
        }
        const double mu = norm > 0.0 ? std::max(cross / norm, 0.0) : 0.0;
        for (size_t j : slots) {
            result.controlPoints[j].mu = mu;
            result.controlPoints[j].aperture.weight = mu;
        }
    }
}

} // namespace

ArcSequenceResult ArcSequencer::conformalSeed(const Stf& sectorStf,
                                              const VmatSettings& settings,
                                              const Machine& machine,
                                              const ArcSequencerOptions& options) {
    const auto sectors = ArcGeometry::sectors(settings);
    const auto cps = ArcGeometry::controlPoints(settings);
    if (sectorStf.getCount() != sectors.size())
        throw std::invalid_argument("ArcSequencer: STF has " + std::to_string(sectorStf.getCount()) +
                                    " beams but the arcs define " + std::to_string(sectors.size()) + " sectors");
    const auto leafBounds = buildMlcLeafBounds(machine.getGeometry());
    if (leafBounds.empty()) throw std::invalid_argument("ArcSequencer: the conformal seed needs an MLC");

    // Every ray of a sector beam (the target projection) gets unit fluence.
    const std::vector<double> weights(sectorStf.getTotalNumOfBixels(), 1.0);
    std::vector<size_t> offset(sectors.size() + 1, 0);
    for (size_t s = 0; s < sectors.size(); ++s)
        offset[s + 1] = offset[s] + sectorStf.getBeam(s)->getTotalNumOfBixels();

    std::vector<std::vector<size_t>> cpsOfSector(sectors.size());
    for (size_t j = 0; j < cps.size(); ++j) cpsOfSector[cps[j].sectorIndex].push_back(j);

    ArcSequenceResult result;
    result.controlPoints.resize(cps.size());
    for (size_t j = 0; j < cps.size(); ++j) result.controlPoints[j].geometry = cps[j];

    std::vector<SectorRows> rows;
    for (size_t s = 0; s < sectors.size(); ++s)
        rows.push_back(resampleToLeafPairs(FluenceMap::fromBeamWeights(*sectorStf.getBeam(s), weights, offset[s]),
                                           leafBounds));
    conformalSequence(result, sectors, cpsOfSector, rows, leafBounds, sectorStf, weights, offset, options);
    computeTiming(result, machine.getConstraints());
    return result;
}

void ArcSequencer::computeTiming(ArcSequenceResult& result, const MachineConstraints& limits) {
    const double vgMax = limits.gantryRotationSpeed[1] > 0.0 ? limits.gantryRotationSpeed[1] : 6.0;
    const double leafMax = limits.leafSpeed[1] > 0.0 ? limits.leafSpeed[1] : 60.0;
    const double rateMax = limits.monitorUnitRate[1] > 0.0 ? limits.monitorUnitRate[1] : 10.0;

    result.totalMU = 0.0;
    result.totalTimeSec = 0.0;
    result.maxLeafSpeed = 0.0;
    result.minGantrySpeed = vgMax;
    result.slowedControlPoints = 0;
    result.meanLeafTravel = 0.0;
    double travelSum = 0.0;
    size_t travelCount = 0;

    const DeliveredControlPoint* prev = nullptr;
    for (auto& cp : result.controlPoints) {
        const double travel = (prev && prev->geometry.arcIndex == cp.geometry.arcIndex) ? maxLeafTravel(*prev, cp) : 0.0;
        if (prev && prev->geometry.arcIndex == cp.geometry.arcIndex &&
            prev->aperture.bankA.size() == cp.aperture.bankA.size())
            for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) {
                if (cp.aperture.bankB[p] <= cp.aperture.bankA[p] && prev->aperture.bankB[p] <= prev->aperture.bankA[p]) continue;
                travelSum += std::abs(cp.aperture.bankA[p] - prev->aperture.bankA[p]) +
                             std::abs(cp.aperture.bankB[p] - prev->aperture.bankB[p]);
                travelCount += 2;
            }
        const double tGantry = cp.geometry.width / vgMax;
        const double tDose = cp.mu / rateMax;
        const double tLeaf = travel / leafMax;
        cp.timeSec = std::max({tGantry, tDose, tLeaf});

        cp.gantrySpeed = cp.timeSec > 0.0 ? cp.geometry.width / cp.timeSec : 0.0;
        cp.doseRate = cp.timeSec > 0.0 ? cp.mu / cp.timeSec : 0.0;
        cp.leafSpeed = cp.timeSec > 0.0 ? travel / cp.timeSec : 0.0;

        if (cp.timeSec > tGantry * (1.0 + 1e-9)) ++result.slowedControlPoints;
        result.totalMU += cp.mu;
        result.totalTimeSec += cp.timeSec;
        result.maxLeafSpeed = std::max(result.maxLeafSpeed, cp.leafSpeed);
        result.minGantrySpeed = std::min(result.minGantrySpeed, cp.gantrySpeed);
        prev = &cp;
    }
    if (travelCount > 0) result.meanLeafTravel = travelSum / static_cast<double>(travelCount);
}

ArcSequenceResult ArcSequencer::resample(const ArcSequenceResult& source,
                                         const std::vector<ArcControlPoint>& geometry,
                                         const MachineConstraints& limits) {
    if (source.controlPoints.empty()) throw std::invalid_argument("ArcSequencer::resample: empty sequence");
    ArcSequenceResult out;
    out.controlPoints.resize(geometry.size());
    for (size_t j = 0; j < geometry.size(); ++j) {
        const ArcControlPoint& target = geometry[j];
        // Source control points of this arc form a contiguous run ordered by offset.
        size_t first = source.controlPoints.size(), last = 0;
        for (size_t k = 0; k < source.controlPoints.size(); ++k) {
            if (source.controlPoints[k].geometry.arcIndex != target.arcIndex) continue;
            first = std::min(first, k);
            last = std::max(last, k);
        }
        if (first > last) throw std::invalid_argument("ArcSequencer::resample: arc missing in the source sequence");

        size_t upper = first;
        while (upper < last && source.controlPoints[upper].geometry.offset < target.offset) ++upper;
        const size_t lower = upper > first ? upper - 1 : upper;
        const auto& a = source.controlPoints[lower];
        const auto& b = source.controlPoints[upper];
        const double span = b.geometry.offset - a.geometry.offset;
        const double t = span > 0.0 ? std::clamp((target.offset - a.geometry.offset) / span, 0.0, 1.0) : 0.0;

        DeliveredControlPoint& cp = out.controlPoints[j];
        cp.geometry = target;
        cp.leafPairBoundariesZ = a.leafPairBoundariesZ;
        cp.aperture.bankA.resize(a.aperture.bankA.size());
        cp.aperture.bankB.resize(a.aperture.bankB.size());
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) {
            cp.aperture.bankA[p] = (1.0 - t) * a.aperture.bankA[p] + t * b.aperture.bankA[p];
            cp.aperture.bankB[p] = (1.0 - t) * a.aperture.bankB[p] + t * b.aperture.bankB[p];
        }
        const double perDegree = (1.0 - t) * a.mu / a.geometry.width + t * b.mu / b.geometry.width;
        cp.mu = perDegree * target.width;
        cp.aperture.weight = cp.mu;
        cp.aperture.segmentIndex = static_cast<int>(j);
    }
    computeTiming(out, limits);
    return out;
}

} // namespace optirad
