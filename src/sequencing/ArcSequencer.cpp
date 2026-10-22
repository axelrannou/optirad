#include "ArcSequencer.hpp"

#include "core/FluenceMap.hpp"
#include "LeafSequencer.hpp"
#include "ApertureFluence.hpp"
#include "MlcLeafBounds.hpp"
#include "SlidingWindow.hpp"

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

// Match physical leaf pairs, including closed pairs.
int matchPair(const DeliveredControlPoint& src, size_t i, const DeliveredControlPoint& dst) {
    if (i + 1 >= src.leafPairBoundariesZ.size() || i >= src.aperture.bankA.size()) return -1;
    const double zc = 0.5 * (src.leafPairBoundariesZ[i] + src.leafPairBoundariesZ[i + 1]);
    const auto& z = dst.leafPairBoundariesZ;
    for (size_t j = 0; j + 1 < z.size() && j < dst.aperture.bankA.size(); ++j) {
        if (zc < z[j] || zc >= z[j + 1]) continue;
        return static_cast<int>(j);
    }
    return -1;
}

// Limits leaf displacement between consecutive control points of every arc to `maxSpeed` at full gantry speed.
// A forward and a backward clamp are averaged (both feasible, so the average is too) so the trajectory
// lags neither the start nor the end of a movement.
void limitLeafSpeed(std::vector<DeliveredControlPoint>& cps, double maxSpeed, double maxGantrySpeed) {
    size_t s = 0;
    while (s < cps.size()) {
        size_t e = s + 1;
        while (e < cps.size() && cps[e].geometry.arcIndex == cps[s].geometry.arcIndex) ++e;

        std::vector<Aperture> fwd(e - s), bwd(e - s);
        for (size_t k = s; k < e; ++k) fwd[k - s] = bwd[k - s] = cps[k].aperture;

        auto clampTo = [](double v, double ref, double lim) { return std::min(std::max(v, ref - lim), ref + lim); };

        for (size_t k = s + 1; k < e; ++k) {
            const double lim = 0.99 * maxSpeed * cps[k].geometry.width / maxGantrySpeed;
            for (size_t i = 0; i < cps[k].aperture.bankA.size(); ++i) {
                const int j = matchPair(cps[k], i, cps[k - 1]);
                if (j < 0) continue;
                fwd[k - s].bankA[i] = clampTo(fwd[k - s].bankA[i], fwd[k - 1 - s].bankA[static_cast<size_t>(j)], lim);
                fwd[k - s].bankB[i] = clampTo(fwd[k - s].bankB[i], fwd[k - 1 - s].bankB[static_cast<size_t>(j)], lim);
            }
        }
        for (size_t k = e - 1; k > s; --k) {
            const double lim = 0.99 * maxSpeed * cps[k].geometry.width / maxGantrySpeed;
            for (size_t i = 0; i < cps[k].aperture.bankA.size(); ++i) {
                const int j = matchPair(cps[k], i, cps[k - 1]);
                if (j < 0) continue;
                const auto jj = static_cast<size_t>(j);
                bwd[k - 1 - s].bankA[jj] = clampTo(bwd[k - 1 - s].bankA[jj], bwd[k - s].bankA[i], lim);
                bwd[k - 1 - s].bankB[jj] = clampTo(bwd[k - 1 - s].bankB[jj], bwd[k - s].bankB[i], lim);
            }
        }
        for (size_t k = s; k < e; ++k) {
            auto& ap = cps[k].aperture;
            for (size_t i = 0; i < ap.bankA.size(); ++i) {
                const double a = 0.5 * (fwd[k - s].bankA[i] + bwd[k - s].bankA[i]);
                const double b = 0.5 * (fwd[k - s].bankB[i] + bwd[k - s].bankB[i]);
                if (b > a) { ap.bankA[i] = a; ap.bankB[i] = b; }
                else { ap.bankA[i] = ap.bankB[i] = 0.5 * (a + b); }
            }
        }
        s = e;
    }
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

    // Ignore negligible fluence so it does not widen the sweep.
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

// Smooths (mass preserving [1 2 1] passes) the rows of a sector's fluence whose upward steps exceed
// `ratio` times the sector's peak fluence, so that no single row sets the MU of the sector.
void limitRowRises(FluenceMap& fluence, double ratio) {
    const int columns = fluence.getNumCols();
    if (ratio <= 0.0 || columns < 3) return;
    const double limit = ratio * fluence.getMaxFluence();
    if (limit <= 0.0) return;
    for (int r = 0; r < fluence.getNumRows(); ++r) {
        std::vector<double> row = fluence.getProfile(r);
        for (int pass = 0; pass < 60 && SlidingWindow::totalRise(row) > limit; ++pass) {
            std::vector<double> next(row.size());
            for (size_t c = 0; c < row.size(); ++c) {
                const double left = row[c == 0 ? 0 : c - 1], right = row[c + 1 == row.size() ? c : c + 1];
                next[c] = 0.25 * left + 0.5 * row[c] + 0.25 * right;
            }
            row = std::move(next);
        }
        for (int c = 0; c < columns; ++c) fluence.setValue(r, c, row[static_cast<size_t>(c)]);
    }
}

// Conformal seed: each leaf pair is open over the extent of the sector fluence above a threshold; the extents are
// interpolated between sector centres so the leaves glide along the arc. MU per control point is the least-squares
// scale that matches the sector's optimized fluence.
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

bool sweepApplicable(const MachineGeometry& mlc, const std::vector<std::vector<size_t>>& cpsOfSector) {
    if (!mlc.interdigitation || buildMlcLeafBounds(mlc).empty()) return false;
    for (const auto& slots : cpsOfSector) if (slots.size() < 2) return false;
    return true;
}

std::vector<double> sweepSequence(ArcSequenceResult& result,
                   const std::vector<ArcSector>& sectors,
                   const std::vector<std::vector<size_t>>& cpsOfSector,
                   const std::vector<SectorRows>& rows,
                   const MachineGeometry& mlc,
                   const MachineConstraints& limits,
                   const ArcSequencerOptions& options) {
    const std::vector<double> leafBounds = buildMlcLeafBounds(mlc);
    const size_t pairs = leafBounds.size() - 1;
    const double leafMax = limits.leafSpeed[1] > 0.0 ? limits.leafSpeed[1] : 60.0;
    const double gantryMax = limits.gantryRotationSpeed[1] > 0.0 ? limits.gantryRotationSpeed[1] : 6.0;
    const double rateMax = limits.monitorUnitRate[1] > 0.0 ? limits.monitorUnitRate[1] : 10.0;
    const double resolution = options.leafPositionResolution;
    const double fraction = std::clamp(options.maxSweepFraction, 0.05, 0.95);
    auto snap = [resolution](double v) { return resolution > 0.0 ? std::round(v / resolution) * resolution : v; };
    const double low = -mlc.jawX1Max, high = mlc.jawX2Max;

    std::vector<double> sectorSpeed(sectors.size(), 0.0);
    size_t arcStart = 0;
    while (arcStart < sectors.size()) {
        size_t arcEnd = arcStart + 1;
        while (arcEnd < sectors.size() && sectors[arcEnd].arcIndex == sectors[arcStart].arcIndex) ++arcEnd;

        for (size_t s = arcStart; s < arcEnd; ++s) {
            const bool forward = ((s - arcStart) % 2) == 0;
            const auto& slots = cpsOfSector[s];
            const int count = static_cast<int>(slots.size());
            const SectorRows& mine = rows[s];
            const SectorRows* before = s > arcStart ? &rows[s - 1] : nullptr;
            const SectorRows* after = s + 1 < arcEnd ? &rows[s + 1] : nullptr;

            // Sweep range per pair in the sweep coordinate y (y = x forward, y = -x backward). Neighbouring
            // sectors share the position where one ends and the next starts.
            std::vector<double> yStart(pairs), yEnd(pairs);
            for (size_t p = 0; p < pairs; ++p) {
                double xStart, xEnd;
                if (forward) {
                    xStart = before ? std::min(before->left[p], mine.left[p]) : mine.left[p];
                    xEnd = after ? std::max(after->right[p], mine.right[p]) : mine.right[p];
                    yStart[p] = xStart; yEnd[p] = xEnd;
                } else {
                    xStart = before ? std::max(before->right[p], mine.right[p]) : mine.right[p];
                    xEnd = after ? std::min(after->left[p], mine.left[p]) : mine.left[p];
                    yStart[p] = -xStart; yEnd[p] = -xEnd;
                }
            }

            std::vector<SlidingWindowRow> sweepRows(pairs);
            std::vector<double> rise(pairs, 0.0), travel(pairs, 0.0);
            double longest = 0.0;
            for (size_t p = 0; p < pairs; ++p) {
                auto& row = sweepRows[p];
                row.pixel = mine.pixel;
                row.start = yStart[p];
                row.end = yEnd[p];
                row.profile = mine.profile[p];
                if (forward) {
                    row.firstCenter = mine.originX;
                } else {
                    std::reverse(row.profile.begin(), row.profile.end());
                    row.firstCenter = -(mine.originX + (mine.columns - 1) * mine.pixel);
                }
                rise[p] = SlidingWindow::totalRise(row.profile);
                travel[p] = std::max(0.0, row.end - row.start);
                longest = std::max(longest, travel[p]);
            }

            // Slice time t (the gantry slows down when needed) fixes the sweep speed v = leafMax * t and the MU per
            // slice m = max over rows of rise / (count - travel / v); the slice time must also allow that MU at the
            // maximum dose rate. The fixed point of t = max(gantry time, m / rate) minimizes the MU.
            const double width = result.controlPoints[slots.front()].geometry.width;
            const double gantryTime = width / gantryMax;
            const double leafTravelPerTime = 0.97 * leafMax;
            const double minTime = std::max(gantryTime, longest / (fraction * count * leafTravelPerTime));
            auto muPerSlotFor = [&](double time) {
                double m = 0.0;
                for (size_t p = 0; p < pairs; ++p)
                    if (rise[p] > 0.0) m = std::max(m, rise[p] / (count - travel[p] / (leafTravelPerTime * time)));
                return m;
            };
            auto excess = [&](double time) { return std::max(gantryTime, muPerSlotFor(time) / rateMax) - time; };
            double sliceTime = minTime;
            if (excess(sliceTime) > 0.0) {
                double high = 2.0 * sliceTime;
                for (int i = 0; i < 60 && excess(high) > 0.0; ++i) high *= 2.0;
                double low = sliceTime;
                for (int i = 0; i < 60; ++i) {
                    const double mid = 0.5 * (low + high);
                    (excess(mid) > 0.0 ? low : high) = mid;
                }
                sliceTime = high;
            }
            const double speed = leafTravelPerTime * sliceTime;
            const double muPerSlot = muPerSlotFor(sliceTime);
            sectorSpeed[s] = speed;

            std::vector<std::vector<double>> bankA(static_cast<size_t>(count), std::vector<double>(pairs));
            auto bankB = bankA;
            std::vector<double> trailing, leading;
            for (size_t p = 0; p < pairs; ++p) {
                SlidingWindow::sampleTips(sweepRows[p], speed, muPerSlot, count, trailing, leading);
                for (int k = 0; k < count; ++k) {
                    const double t = trailing[static_cast<size_t>(k)], l = leading[static_cast<size_t>(k)];
                    // Forward: the leading leaf is the right bank. Backward the axis is mirrored.
                    const double a = forward ? t : -l, b = forward ? l : -t;
                    bankA[static_cast<size_t>(k)][p] = snap(std::clamp(a, low, high));
                    bankB[static_cast<size_t>(k)][p] = snap(std::clamp(b, low, high));
                }
            }

            for (int k = 0; k < count; ++k) {
                DeliveredControlPoint& cp = result.controlPoints[slots[static_cast<size_t>(k)]];
                cp.leafPairBoundariesZ = leafBounds;
                cp.aperture.bankA = std::move(bankA[static_cast<size_t>(k)]);
                cp.aperture.bankB = std::move(bankB[static_cast<size_t>(k)]);
                cp.mu = muPerSlot;
                cp.aperture.weight = muPerSlot;
                cp.aperture.segmentIndex = static_cast<int>(slots[static_cast<size_t>(k)]);
            }
        }
        arcStart = arcEnd;
    }
    return sectorSpeed;
}

// Refines the leaf tips of every control point so that the fluence delivered by the moving leaves (summed over
// the control points of a sector, on the sector's optimization beam) matches the optimized fluence.
// Projected gradient with Barzilai-Borwein steps; tips stay within the jaws, keep left <= right and move at
// most the sector's sweep speed per slice.
void fitToFluence(ArcSequenceResult& result,
                  const std::vector<std::vector<size_t>>& cpsOfSector,
                  const std::vector<double>& sectorSpeed,
                  const Stf& fmoStf,
                  const std::vector<double>& weights,
                  const std::vector<size_t>& offset,
                  const MachineGeometry& mlc,
                  int iterations) {
    const size_t count = result.controlPoints.size();
    if (count == 0 || iterations <= 0) return;
    const size_t pairs = result.controlPoints.front().aperture.bankA.size();
    const double low = -mlc.jawX1Max, high = mlc.jawX2Max;

    std::vector<Aperture> tips(count), gradient(count), candidate(count), previousTips, previousGradient;
    for (size_t j = 0; j < count; ++j) {
        tips[j] = result.controlPoints[j].aperture;
        gradient[j] = tips[j];
    }
    auto sameArc = [&](size_t a, size_t b) {
        return result.controlPoints[a].geometry.arcIndex == result.controlPoints[b].geometry.arcIndex;
    };
    auto neighbour = [&](const std::vector<Aperture>& set, size_t j, int step) -> const Aperture* {
        if ((step < 0 && j == 0) || (step > 0 && j + 1 >= count)) return nullptr;
        const size_t other = step < 0 ? j - 1 : j + 1;
        return sameArc(j, other) ? &set[other] : nullptr;
    };
    auto neighbourIndex = [&](size_t j, int step) -> long {
        if ((step < 0 && j == 0) || (step > 0 && j + 1 >= count)) return -1;
        const size_t other = step < 0 ? j - 1 : j + 1;
        return sameArc(j, other) ? static_cast<long>(other) : -1;
    };

    auto evaluate = [&](const std::vector<Aperture>& set, std::vector<Aperture>* grad) {
        std::vector<std::vector<double>> residual(fmoStf.getCount());
        for (size_t s = 0; s < fmoStf.getCount(); ++s) {
            const size_t rays = fmoStf.getBeam(s)->getNumOfRays();
            residual[s].resize(rays);
            for (size_t r = 0; r < rays; ++r) residual[s][r] = -weights[offset[s] + r];
        }
        for (size_t s = 0; s < cpsOfSector.size(); ++s) {
            for (size_t j : cpsOfSector[s]) {
                const double mu = result.controlPoints[j].mu;
                const auto open = ApertureFluence::integratedOpenFractions(*fmoStf.getBeam(s),
                    neighbour(set, j, -1), set[j], neighbour(set, j, 1), result.controlPoints[j].leafPairBoundariesZ);
                for (size_t r = 0; r < open.size(); ++r) residual[s][r] += mu * open[r];
            }
        }
        double value = 0.0;
        for (const auto& sector : residual) for (double v : sector) value += 0.5 * v * v;
        if (!grad) return value;
        for (auto& g : *grad) {
            std::fill(g.bankA.begin(), g.bankA.end(), 0.0);
            std::fill(g.bankB.begin(), g.bankB.end(), 0.0);
        }
        // Control points j, j+3, ... write to disjoint apertures, so each residue class can run in parallel.
        for (size_t color = 0; color < 3; ++color) {
            #pragma omp parallel for schedule(dynamic)
            for (long j = static_cast<long>(color); j < static_cast<long>(count); j += 3) {
                const size_t index = static_cast<size_t>(j);
                const size_t s = result.controlPoints[index].geometry.sectorIndex;
                const double mu = result.controlPoints[index].mu;
                std::vector<double> rayWeight(residual[s].size());
                for (size_t r = 0; r < rayWeight.size(); ++r) rayWeight[r] = mu * residual[s][r];
                const long before = neighbourIndex(index, -1), after = neighbourIndex(index, 1);
                ApertureFluence::accumulateIntegratedTipGradient(*fmoStf.getBeam(s),
                    before >= 0 ? &set[static_cast<size_t>(before)] : nullptr, set[index],
                    after >= 0 ? &set[static_cast<size_t>(after)] : nullptr,
                    result.controlPoints[index].leafPairBoundariesZ, rayWeight,
                    before >= 0 ? &(*grad)[static_cast<size_t>(before)] : nullptr, (*grad)[index],
                    after >= 0 ? &(*grad)[static_cast<size_t>(after)] : nullptr);
            }
        }
        return value;
    };

    auto project = [&](std::vector<Aperture>& set) {
        for (size_t j = 0; j < count; ++j) {
            const double speed = sectorSpeed[result.controlPoints[j].geometry.sectorIndex];
            const bool hasPrevious = neighbourIndex(j, -1) >= 0;
            for (size_t p = 0; p < pairs; ++p) {
                double a = std::clamp(set[j].bankA[p], low, high), b = std::clamp(set[j].bankB[p], low, high);
                if (hasPrevious) {
                    const Aperture& prev = set[j - 1];
                    a = std::clamp(a, prev.bankA[p] - speed, prev.bankA[p] + speed);
                    b = std::clamp(b, prev.bankB[p] - speed, prev.bankB[p] + speed);
                }
                if (a > b) a = b = 0.5 * (a + b);
                set[j].bankA[p] = a;
                set[j].bankB[p] = b;
            }
        }
    };

    // Delivered fluence integral per sector should match the target before and after the fit.
    auto rescaleMu = [&]() {
        for (size_t s = 0; s < cpsOfSector.size(); ++s) {
            std::vector<double> delivered(fmoStf.getBeam(s)->getNumOfRays(), 0.0);
            for (size_t j : cpsOfSector[s]) {
                const auto open = ApertureFluence::integratedOpenFractions(*fmoStf.getBeam(s),
                    neighbour(tips, j, -1), tips[j], neighbour(tips, j, 1), result.controlPoints[j].leafPairBoundariesZ);
                for (size_t r = 0; r < open.size(); ++r) delivered[r] += result.controlPoints[j].mu * open[r];
            }
            double cross = 0.0, norm = 0.0;
            for (size_t r = 0; r < delivered.size(); ++r) {
                cross += delivered[r] * weights[offset[s] + r];
                norm += delivered[r] * delivered[r];
            }
            if (norm <= 0.0) continue;
            const double scale = std::clamp(cross / norm, 0.5, 2.0);
            for (size_t j : cpsOfSector[s]) result.controlPoints[j].mu *= scale;
        }
    };

    project(tips);
    rescaleMu();
    double value = evaluate(tips, &gradient);
    double step = 0.0;
    for (int iteration = 0; iteration < iterations; ++iteration) {
        double maxGradient = 0.0;
        for (const auto& g : gradient) {
            for (double v : g.bankA) maxGradient = std::max(maxGradient, std::abs(v));
            for (double v : g.bankB) maxGradient = std::max(maxGradient, std::abs(v));
        }
        if (maxGradient <= 0.0) break;
        if (iteration == 0 || step <= 0.0) step = 1.0 / maxGradient;   // first move: 1 mm

        bool accepted = false;
        double trial = step;
        for (int attempt = 0; attempt < 12; ++attempt, trial *= 0.5) {
            for (size_t j = 0; j < count; ++j) {
                candidate[j] = tips[j];
                for (size_t p = 0; p < pairs; ++p) {
                    candidate[j].bankA[p] -= trial * gradient[j].bankA[p];
                    candidate[j].bankB[p] -= trial * gradient[j].bankB[p];
                }
            }
            project(candidate);
            const double next = evaluate(candidate, nullptr);
            if (!std::isfinite(next) || next >= value) continue;
            previousTips = std::move(tips);
            previousGradient = std::move(gradient);
            tips = candidate;
            gradient = tips;
            value = evaluate(tips, &gradient);
            accepted = true;
            break;
        }
        if (!accepted) break;

        double ss = 0.0, sy = 0.0;
        for (size_t j = 0; j < count; ++j)
            for (size_t p = 0; p < pairs; ++p) {
                const double da = tips[j].bankA[p] - previousTips[j].bankA[p];
                const double db = tips[j].bankB[p] - previousTips[j].bankB[p];
                ss += da * da + db * db;
                sy += da * (gradient[j].bankA[p] - previousGradient[j].bankA[p]) +
                      db * (gradient[j].bankB[p] - previousGradient[j].bankB[p]);
            }
        step = (sy > 1e-30 && ss > 0.0) ? std::clamp(ss / sy, 1e-3 / maxGradient, 50.0 / maxGradient) : step;
    }
    rescaleMu();
    for (size_t j = 0; j < count; ++j) {
        result.controlPoints[j].aperture.bankA = tips[j].bankA;
        result.controlPoints[j].aperture.bankB = tips[j].bankB;
        result.controlPoints[j].aperture.weight = result.controlPoints[j].mu;
    }
}

} // namespace

ArcSequenceResult ArcSequencer::sequence(const std::vector<double>& weights,
                                         const Stf& fmoStf,
                                         const VmatSettings& settings,
                                         const Machine& machine,
                                         const ArcSequencerOptions& options) {
    const auto sectors = ArcGeometry::sectors(settings);
    const auto cps = ArcGeometry::controlPoints(settings);
    if (fmoStf.getCount() != sectors.size())
        throw std::invalid_argument("ArcSequencer: STF has " + std::to_string(fmoStf.getCount()) +
                                    " beams but the arcs define " + std::to_string(sectors.size()) + " sectors");

    std::vector<size_t> offset(sectors.size() + 1, 0);
    for (size_t s = 0; s < sectors.size(); ++s)
        offset[s + 1] = offset[s] + fmoStf.getBeam(s)->getTotalNumOfBixels();

    std::vector<std::vector<size_t>> cpsOfSector(sectors.size());
    for (size_t j = 0; j < cps.size(); ++j) cpsOfSector[cps[j].sectorIndex].push_back(j);

    ArcSequenceResult result;
    result.controlPoints.resize(cps.size());
    for (size_t j = 0; j < cps.size(); ++j) result.controlPoints[j].geometry = cps[j];

    const MachineGeometry& mlc = machine.getGeometry();

    if (options.method == ArcSequencerOptions::Method::Conformal) {
        const auto leafBounds = buildMlcLeafBounds(mlc);
        if (leafBounds.empty()) throw std::invalid_argument("ArcSequencer: conformal sequencing needs an MLC");
        std::vector<SectorRows> rows;
        for (size_t s = 0; s < sectors.size(); ++s)
            rows.push_back(resampleToLeafPairs(FluenceMap::fromBeamWeights(*fmoStf.getBeam(s), weights, offset[s]),
                                               leafBounds));
        conformalSequence(result, sectors, cpsOfSector, rows, leafBounds, fmoStf, weights, offset, options);
        computeTiming(result, machine.getConstraints());
        return result;
    }

    const bool canSweep = sweepApplicable(mlc, cpsOfSector);
    if (options.method == ArcSequencerOptions::Method::Sweep && !canSweep)
        throw std::invalid_argument("ArcSequencer: sweep sequencing needs an interdigitating MLC and at least "
                                    "two control points per sector");
    if (options.method != ArcSequencerOptions::Method::LevelSet && canSweep) {
        const auto leafBounds = buildMlcLeafBounds(mlc);
        std::vector<SectorRows> rows;
        std::vector<double> target = weights;   // the fluence the sweep delivers (rows limited in modulation)
        for (size_t s = 0; s < sectors.size(); ++s) {
            FluenceMap fluence = FluenceMap::fromBeamWeights(*fmoStf.getBeam(s), weights, offset[s]);
            limitRowRises(fluence, options.maxRiseRatio);
            fluence.toBeamWeights(*fmoStf.getBeam(s), target, offset[s]);
            rows.push_back(resampleToLeafPairs(fluence, leafBounds));
        }
        const auto sectorSpeed = sweepSequence(result, sectors, cpsOfSector, rows, mlc, machine.getConstraints(),
                                               options);
        fitToFluence(result, cpsOfSector, sectorSpeed, fmoStf, target, offset, mlc, options.fitIterations);
        const double resolution = options.leafPositionResolution;
        if (resolution > 0.0)
            for (auto& cp : result.controlPoints) {
                for (double& v : cp.aperture.bankA) v = std::round(v / resolution) * resolution;
                for (double& v : cp.aperture.bankB) v = std::round(v / resolution) * resolution;
            }
        computeTiming(result, machine.getConstraints());
        return result;
    }

    size_t currentArc = static_cast<size_t>(-1);
    size_t sectorInArc = 0;
    for (size_t s = 0; s < sectors.size(); ++s) {
        if (sectors[s].arcIndex != currentArc) { currentArc = sectors[s].arcIndex; sectorInArc = 0; }
        const bool reverseOrder = (sectorInArc++ % 2) == 1;

        const auto& slots = cpsOfSector[s];
        if (slots.empty()) continue;
        const int nCP = static_cast<int>(slots.size());

        const Beam* beam = fmoStf.getBeam(s);
        const FluenceMap fluence = FluenceMap::fromBeamWeights(*beam, weights, offset[s]);

        LeafSequencerOptions lso;
        lso.numLevels = std::max(1, std::min(options.numLevels, nCP));
        lso.minSegmentMU = options.minSegmentMU;
        lso.leafPositionResolution = options.leafPositionResolution;
        LeafSequenceResult seq = LeafSequencer::sequenceBeam(fluence, mlc, lso);

        std::vector<double> zBounds = seq.leafPairBoundariesZ;
        if (zBounds.empty() && fluence.getNumRows() > 0) {
            for (int r = 0; r <= fluence.getNumRows(); ++r)
                zBounds.push_back(fluence.getOriginZ() - 0.5 * fluence.getBixelWidth() + r * fluence.getBixelWidth());
        }

        std::vector<Aperture> segments = seq.segments;
        if (reverseOrder) std::reverse(segments.begin(), segments.end());
        const int S = static_cast<int>(segments.size());

        std::vector<int> perSegment(static_cast<size_t>(std::max(S, 1)), 0);
        for (int k = 0; k < nCP && S > 0; ++k) ++perSegment[static_cast<size_t>(k * S / nCP)];

        for (int k = 0; k < nCP; ++k) {
            DeliveredControlPoint& cp = result.controlPoints[slots[static_cast<size_t>(k)]];
            cp.leafPairBoundariesZ = zBounds;
            if (S == 0) {
                const size_t pairs = zBounds.empty() ? 0 : zBounds.size() - 1;
                cp.aperture.bankA.assign(pairs, 0.0);
                cp.aperture.bankB.assign(pairs, 0.0);
                cp.mu = 0.0;
            } else {
                const int idx = k * S / nCP;
                cp.aperture = segments[static_cast<size_t>(idx)];
                cp.mu = segments[static_cast<size_t>(idx)].weight / perSegment[static_cast<size_t>(idx)];
            }
            cp.aperture.weight = cp.mu;
            cp.aperture.segmentIndex = static_cast<int>(slots[static_cast<size_t>(k)]);
        }
    }

    std::vector<double> commonBounds = buildMlcLeafBounds(mlc);
    if (commonBounds.empty()) {
        for (const auto& cp : result.controlPoints)
            commonBounds.insert(commonBounds.end(), cp.leafPairBoundariesZ.begin(), cp.leafPairBoundariesZ.end());
    }
    std::sort(commonBounds.begin(), commonBounds.end());
    commonBounds.erase(std::unique(commonBounds.begin(), commonBounds.end()), commonBounds.end());
    for (auto& cp : result.controlPoints) {
        Aperture mapped = cp.aperture;
        mapped.bankA.assign(commonBounds.empty() ? 0 : commonBounds.size() - 1, 0.0);
        mapped.bankB = mapped.bankA;
        for (size_t pair = 0; pair < mapped.bankA.size(); ++pair) {
            const double center = 0.5 * (commonBounds[pair] + commonBounds[pair + 1]);
            for (size_t old = 0; old < cp.aperture.bankA.size(); ++old) {
                if (center < cp.leafPairBoundariesZ[old] || center >= cp.leafPairBoundariesZ[old + 1]) continue;
                mapped.bankA[pair] = cp.aperture.bankA[old];
                mapped.bankB[pair] = cp.aperture.bankB[old];
                break;
            }
        }
        cp.aperture = std::move(mapped);
        cp.leafPairBoundariesZ = commonBounds;
    }
    const MachineConstraints& limits = machine.getConstraints();
    if (options.enforceLeafSpeed && limits.leafSpeed[1] > 0.0 && limits.gantryRotationSpeed[1] > 0.0)
        limitLeafSpeed(result.controlPoints, limits.leafSpeed[1], limits.gantryRotationSpeed[1]);

    computeTiming(result, limits);
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
