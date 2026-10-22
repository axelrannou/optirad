#include "ApertureFluence.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace optirad {

std::vector<double> ApertureFluence::bixelOpenFractions(const Beam& beam,
                                                        const Aperture& aperture,
                                                        const std::vector<double>& zBounds) {
    const size_t numRays = beam.getNumOfRays();
    std::vector<double> open(numRays, 0.0);

    const size_t pairs = std::min(aperture.bankA.size(), aperture.bankB.size());
    if (pairs == 0 || zBounds.size() < pairs + 1) return open;

    const double bw = beam.getBixelWidth();
    const double half = 0.5 * bw;

    for (size_t r = 0; r < numRays; ++r) {
        const auto& pos = beam.getRay(r)->getRayPosBev();
        const double x0 = pos[0] - half, x1 = pos[0] + half;
        const double z0 = pos[2] - half, z1 = pos[2] + half;

        double fraction = 0.0;
        for (size_t lp = 0; lp < pairs; ++lp) {
            const double overlapZ = std::min(z1, zBounds[lp + 1]) - std::max(z0, zBounds[lp]);
            if (overlapZ <= 0.0) continue;
            const double overlapX = std::min(x1, aperture.bankB[lp]) - std::max(x0, aperture.bankA[lp]);
            if (overlapX <= 0.0) continue;
            fraction += (overlapZ / bw) * (overlapX / bw);
        }
        open[r] = std::min(1.0, fraction);
    }
    return open;
}

namespace {

// Time integral over t in [0,1] of the overlap between the pixel [x0,x1] and the window [A(t),B(t)], with
// A and B linear between (a0,a1) and (b0,b1). `grad` receives the derivatives with respect to a0,a1,b0,b1.
double integrateOverlap(double x0, double x1, double a0, double a1, double b0, double b1,
                        std::array<double, 4>* grad) {
    std::array<double, 6> knots = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    size_t count = 2;
    auto addCrossing = [&](double from, double to, double level) {
        if (from == to) return;
        const double t = (level - from) / (to - from);
        if (t > 0.0 && t < 1.0) knots[count++] = t;
    };
    addCrossing(a0, a1, x0);
    addCrossing(a0, a1, x1);
    addCrossing(b0, b1, x0);
    addCrossing(b0, b1, x1);
    std::sort(knots.begin(), knots.begin() + count);

    auto lerp = [](double from, double to, double t) { return from + (to - from) * t; };
    auto gap = [&](double t) {
        return std::min(x1, lerp(b0, b1, t)) - std::max(x0, lerp(a0, a1, t));
    };

    double total = 0.0;
    if (grad) grad->fill(0.0);
    for (size_t i = 0; i + 1 < count; ++i) {
        const double t0 = knots[i], t1 = knots[i + 1], width = t1 - t0;
        if (width <= 0.0) continue;
        const double g0 = gap(t0), g1 = gap(t1);
        if (g0 <= 0.0 && g1 <= 0.0) continue;
        double lo = t0, hi = t1;
        if (g0 > 0.0 && g1 > 0.0) total += 0.5 * (g0 + g1) * width;
        else if (g0 > 0.0) { hi = t0 + width * g0 / (g0 - g1); total += 0.5 * g0 * (hi - t0); }
        else { lo = t1 - width * g1 / (g1 - g0); total += 0.5 * g1 * (t1 - lo); }
        if (!grad) continue;
        const double tm = 0.5 * (lo + hi);
        const double am = lerp(a0, a1, tm), bm = lerp(b0, b1, tm);
        const double linear = 0.5 * (hi * hi - lo * lo);
        const double early = (hi - lo) - linear;   // integral of (1 - t)
        if (am > x0 && am < x1) { (*grad)[0] -= early; (*grad)[1] -= linear; }
        if (bm > x0 && bm < x1) { (*grad)[2] += early; (*grad)[3] += linear; }
    }
    return total;
}

} // namespace

std::vector<double> ApertureFluence::integratedOpenFractions(const Beam& beam,
                                                             const Aperture* previous,
                                                             const Aperture& current,
                                                             const Aperture* next,
                                                             const std::vector<double>& zBounds) {
    const size_t numRays = beam.getNumOfRays();
    std::vector<double> open(numRays, 0.0);
    const size_t pairs = std::min(current.bankA.size(), current.bankB.size());
    if (pairs == 0 || zBounds.size() < pairs + 1) return open;
    const Aperture& before = previous ? *previous : current;
    const Aperture& after = next ? *next : current;
    if (before.bankA.size() < pairs || before.bankB.size() < pairs ||
        after.bankA.size() < pairs || after.bankB.size() < pairs)
        return bixelOpenFractions(beam, current, zBounds);

    const double bw = beam.getBixelWidth();
    const double half = 0.5 * bw;
    for (size_t r = 0; r < numRays; ++r) {
        const auto& pos = beam.getRay(r)->getRayPosBev();
        const double x0 = pos[0] - half, x1 = pos[0] + half;
        const double z0 = pos[2] - half, z1 = pos[2] + half;
        size_t pair = static_cast<size_t>(std::upper_bound(zBounds.begin(), zBounds.begin() + pairs + 1, z0) -
                                          zBounds.begin());
        pair = pair > 0 ? pair - 1 : 0;
        double fraction = 0.0;
        for (; pair < pairs && zBounds[pair] < z1; ++pair) {
            const double overlapZ = std::min(z1, zBounds[pair + 1]) - std::max(z0, zBounds[pair]);
            if (overlapZ <= 0.0) continue;
            const double a = current.bankA[pair], b = current.bankB[pair];
            const double aBefore = 0.5 * (before.bankA[pair] + a), bBefore = 0.5 * (before.bankB[pair] + b);
            const double aAfter = 0.5 * (after.bankA[pair] + a), bAfter = 0.5 * (after.bankB[pair] + b);
            const double motion = 0.5 * (integrateOverlap(x0, x1, aBefore, a, bBefore, b, nullptr) +
                                         integrateOverlap(x0, x1, a, aAfter, b, bAfter, nullptr));
            fraction += (overlapZ / bw) * (motion / bw);
        }
        open[r] = std::min(1.0, fraction);
    }
    return open;
}

void ApertureFluence::accumulateIntegratedTipGradient(const Beam& beam,
                                                      const Aperture* previous,
                                                      const Aperture& current,
                                                      const Aperture* next,
                                                      const std::vector<double>& zBounds,
                                                      const std::vector<double>& rayWeight,
                                                      Aperture* gradPrevious,
                                                      Aperture& gradCurrent,
                                                      Aperture* gradNext) {
    const size_t pairs = std::min(current.bankA.size(), current.bankB.size());
    if (pairs == 0 || zBounds.size() < pairs + 1) return;
    const Aperture& before = previous ? *previous : current;
    const Aperture& after = next ? *next : current;
    Aperture& gBefore = (previous && gradPrevious) ? *gradPrevious : gradCurrent;
    Aperture& gAfter = (next && gradNext) ? *gradNext : gradCurrent;

    const double bw = beam.getBixelWidth();
    const double half = 0.5 * bw;
    std::array<double, 4> gFirst{}, gSecond{};
    for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
        const double weight = rayWeight[r];
        if (weight == 0.0) continue;
        const auto& pos = beam.getRay(r)->getRayPosBev();
        const double x0 = pos[0] - half, x1 = pos[0] + half;
        const double z0 = pos[2] - half, z1 = pos[2] + half;
        size_t pair = static_cast<size_t>(std::upper_bound(zBounds.begin(), zBounds.begin() + pairs + 1, z0) -
                                          zBounds.begin());
        pair = pair > 0 ? pair - 1 : 0;
        for (; pair < pairs && zBounds[pair] < z1; ++pair) {
            const double overlapZ = std::min(z1, zBounds[pair + 1]) - std::max(z0, zBounds[pair]);
            if (overlapZ <= 0.0) continue;
            const double scale = weight * 0.5 * overlapZ / (bw * bw);
            const double a = current.bankA[pair], b = current.bankB[pair];
            const double aBefore = 0.5 * (before.bankA[pair] + a), bBefore = 0.5 * (before.bankB[pair] + b);
            const double aAfter = 0.5 * (after.bankA[pair] + a), bAfter = 0.5 * (after.bankB[pair] + b);
            integrateOverlap(x0, x1, aBefore, a, bBefore, b, &gFirst);
            integrateOverlap(x0, x1, a, aAfter, b, bAfter, &gSecond);
            // First half: start = (before + current) / 2, end = current.
            gBefore.bankA[pair] += scale * 0.5 * gFirst[0];
            gradCurrent.bankA[pair] += scale * (0.5 * gFirst[0] + gFirst[1]);
            gBefore.bankB[pair] += scale * 0.5 * gFirst[2];
            gradCurrent.bankB[pair] += scale * (0.5 * gFirst[2] + gFirst[3]);
            // Second half: start = current, end = (current + after) / 2.
            gradCurrent.bankA[pair] += scale * (gSecond[0] + 0.5 * gSecond[1]);
            gAfter.bankA[pair] += scale * 0.5 * gSecond[1];
            gradCurrent.bankB[pair] += scale * (gSecond[2] + 0.5 * gSecond[3]);
            gAfter.bankB[pair] += scale * 0.5 * gSecond[3];
        }
    }
}

} // namespace optirad
