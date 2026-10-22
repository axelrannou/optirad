#include <gtest/gtest.h>

#include "core/ArcPlan.hpp"
#include "core/Beam.hpp"
#include "core/FluenceMap.hpp"
#include "core/Ray.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "sequencing/ApertureFluence.hpp"
#include "sequencing/ArcSequencer.hpp"
#include "sequencing/SlidingWindow.hpp"
#include "validation/WaterPhantom.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

namespace optirad::tests {

// ============================================================================
// Arc geometry
// ============================================================================

TEST(ArcGeometryTest, SpanOfClockwiseAndCounterClockwiseArcs) {
    ArcSpec cw;  cw.startAngle = 181; cw.stopAngle = 179; cw.clockwise = true;
    ArcSpec ccw; ccw.startAngle = 179; ccw.stopAngle = 181; ccw.clockwise = false;
    EXPECT_NEAR(cw.span(), 358.0, 1e-9);
    EXPECT_NEAR(ccw.span(), 358.0, 1e-9);

    ArcSpec quarter; quarter.startAngle = 350; quarter.stopAngle = 80; quarter.clockwise = true;
    EXPECT_NEAR(quarter.span(), 90.0, 1e-9);  // wraps through 0
    EXPECT_NEAR(quarter.angleAt(20.0), 10.0, 1e-9);

    ArcSpec full; full.startAngle = 0; full.stopAngle = 0;
    EXPECT_NEAR(full.span(), 360.0, 1e-9);
}

TEST(ArcGeometryTest, SectorsAndControlPointsCoverTheArcWithoutGaps) {
    VmatSettings s;
    s.arcs[0].startAngle = 0; s.arcs[0].stopAngle = 120; s.arcs[0].clockwise = true;
    s.optimizationSpacingDeg = 40.0;
    s.controlPointSpacingDeg = 10.0;

    const auto sectors = ArcGeometry::sectors(s);
    ASSERT_EQ(sectors.size(), 3u);
    EXPECT_NEAR(sectors[0].centerAngle, 20.0, 1e-9);
    EXPECT_NEAR(sectors[2].centerAngle, 100.0, 1e-9);

    const auto cps = ArcGeometry::controlPoints(s);
    ASSERT_EQ(cps.size(), 12u);
    double total = 0.0;
    for (const auto& cp : cps) total += cp.width;
    EXPECT_NEAR(total, 120.0, 1e-9);
    EXPECT_EQ(cps.front().sectorIndex, 0u);
    EXPECT_EQ(cps.back().sectorIndex, 2u);
    for (size_t i = 1; i < cps.size(); ++i) EXPECT_GE(cps[i].sectorIndex, cps[i - 1].sectorIndex);
}

TEST(ArcGeometryTest, TwoArcsKeepTheirOwnSectorsAndCouchAngles) {
    VmatSettings s;
    s.arcs = {ArcSpec{}, ArcSpec{}};
    s.arcs[0].startAngle = 0;   s.arcs[0].stopAngle = 90;  s.arcs[0].clockwise = true;
    s.arcs[1].startAngle = 90;  s.arcs[1].stopAngle = 0;   s.arcs[1].clockwise = false;
    s.arcs[1].couchAngle = 10.0;
    s.optimizationSpacingDeg = 30.0;
    s.controlPointSpacingDeg = 10.0;

    const auto sectors = ArcGeometry::sectors(s);
    ASSERT_EQ(sectors.size(), 6u);
    EXPECT_EQ(sectors[2].arcIndex, 0u);
    EXPECT_EQ(sectors[3].arcIndex, 1u);
    EXPECT_NEAR(sectors[3].centerAngle, 75.0, 1e-9);  // second arc runs backwards
    EXPECT_NEAR(sectors[3].couchAngle, 10.0, 1e-9);

    const auto cps = ArcGeometry::controlPoints(s);
    ASSERT_EQ(cps.size(), 18u);
    EXPECT_EQ(cps[8].arcIndex, 0u);
    EXPECT_EQ(cps[9].arcIndex, 1u);
    EXPECT_EQ(cps[9].sectorIndex, 3u);
}

// ============================================================================
// Aperture -> bixel fractions
// ============================================================================

namespace {

Beam makeBeam(int cols, int rows, double bw) {
    Beam beam;
    beam.setBixelWidth(bw);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) {
            Ray ray;
            ray.setRayPosBev(Vec3{(c - (cols - 1) / 2.0) * bw, 0.0, (r - (rows - 1) / 2.0) * bw});
            beam.addRay(ray);
        }
    return beam;
}

} // namespace

TEST(ApertureFluenceTest, OpenClosedAndPartialBixels) {
    const double bw = 10.0;
    const Beam beam = makeBeam(3, 1, bw);              // bixels centred at x = -10, 0, 10; z = 0
    std::vector<double> z = {-5.0, 5.0};               // one leaf pair covering the row

    Aperture a;
    a.bankA = {-5.0};   // opens x in [-5, 10): bixel 0 is half open (x in [-15,-5] overlaps [-5,-5] = 0)
    a.bankB = {10.0};
    auto f = ApertureFluence::bixelOpenFractions(beam, a, z);
    ASSERT_EQ(f.size(), 3u);
    EXPECT_NEAR(f[0], 0.0, 1e-12);
    EXPECT_NEAR(f[1], 1.0, 1e-12);   // x in [-5, 5] fully open
    EXPECT_NEAR(f[2], 0.5, 1e-12);   // x in [5, 15] open up to 10

    Aperture closed;
    closed.bankA = {0.0};
    closed.bankB = {0.0};
    for (double v : ApertureFluence::bixelOpenFractions(beam, closed, z)) EXPECT_EQ(v, 0.0);
}

TEST(ApertureFluenceTest, BixelSpanningTwoLeafPairsAveragesThem) {
    const double bw = 10.0;
    const Beam beam = makeBeam(1, 1, bw);               // one bixel at the origin, z in [-5, 5]
    std::vector<double> z = {-5.0, 0.0, 5.0};           // two 5 mm leaf pairs

    Aperture a;
    a.bankA = {-5.0, 0.0};   // lower pair open across the bixel, upper pair closed
    a.bankB = {5.0, 0.0};
    auto f = ApertureFluence::bixelOpenFractions(beam, a, z);
    EXPECT_NEAR(f[0], 0.5, 1e-12);
}


TEST(ApertureFluenceTest, IntegratedMotionAveragesTheLeafPathOverTheSlice) {
    const Beam beam = makeBeam(1, 1, 10.0);             // one bixel, x in [-5, 5]
    std::vector<double> z = {-5.0, 5.0};
    auto at = [](double a, double b) { Aperture ap; ap.bankA = {a}; ap.bankB = {b}; return ap; };

    // The right tip follows -7.5 -> 0 -> 7.5 over the slice (previous -15, current 0, next 15).
    const Aperture previous = at(-20.0, -15.0), current = at(-20.0, 0.0), next = at(-20.0, 15.0);
    const double open = ApertureFluence::integratedOpenFractions(beam, &previous, current, &next, z)[0];
    double sampled = 0.0;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        const double u = (i + 0.5) / n;                          // slice position in [0, 1]
        const double tip = u < 0.5 ? -7.5 + 15.0 * u : 15.0 * (u - 0.5);
        sampled += std::max(0.0, std::min(5.0, tip) - (-5.0)) / 10.0 / n;
    }
    EXPECT_NEAR(open, sampled, 1e-4);

    // Without neighbours the slice is static.
    const Aperture still = at(-5.0, 2.0);
    EXPECT_NEAR(ApertureFluence::integratedOpenFractions(beam, nullptr, still, nullptr, z)[0],
                ApertureFluence::bixelOpenFractions(beam, still, z)[0], 1e-12);
}

TEST(ApertureFluenceTest, IntegratedTipGradientMatchesFiniteDifferences) {
    const Beam beam = makeBeam(5, 2, 4.0);
    std::vector<double> z = {-4.0, 0.0, 4.0};
    auto make = [](double a0, double b0, double a1, double b1) {
        Aperture ap; ap.bankA = {a0, a1}; ap.bankB = {b0, b1}; return ap; };
    Aperture previous = make(-7.3, -1.1, -3.9, 2.4), current = make(-4.2, 3.3, -6.1, 5.2), next = make(-0.7, 7.7, -2.2, 6.1);
    std::vector<double> weight(beam.getNumOfRays());
    for (size_t r = 0; r < weight.size(); ++r) weight[r] = 1.0 + 0.37 * static_cast<double>(r % 4);
    auto objective = [&](const Aperture& p, const Aperture& c, const Aperture& n) {
        const auto open = ApertureFluence::integratedOpenFractions(beam, &p, c, &n, z);
        double value = 0.0;
        for (size_t r = 0; r < open.size(); ++r) value += weight[r] * open[r];
        return value;
    };
    Aperture gp, gc, gn;
    for (auto* g : {&gp, &gc, &gn}) { g->bankA.assign(2, 0.0); g->bankB.assign(2, 0.0); }
    ApertureFluence::accumulateIntegratedTipGradient(beam, &previous, current, &next, z, weight, &gp, gc, &gn);

    const double h = 1e-5;
    for (size_t pair = 0; pair < 2; ++pair) {
        for (int which = 0; which < 3; ++which) {
            for (int bank = 0; bank < 2; ++bank) {
                auto shifted = [&](double delta) {
                    Aperture p = previous, c = current, n = next;
                    Aperture& target = which == 0 ? p : (which == 1 ? c : n);
                    (bank == 0 ? target.bankA : target.bankB)[pair] += delta;
                    return objective(p, c, n);
                };
                const double numeric = (shifted(h) - shifted(-h)) / (2.0 * h);
                const Aperture& g = which == 0 ? gp : (which == 1 ? gc : gn);
                const double analytic = (bank == 0 ? g.bankA : g.bankB)[pair];
                EXPECT_NEAR(analytic, numeric, 1e-6) << "pair " << pair << " cp " << which << " bank " << bank;
            }
        }
    }
}

TEST(SlidingWindowTest, FineSlotsReproduceAnArbitraryProfile) {
    const std::vector<double> profile = {0.0, 2.0, 5.0, 1.0, 0.0, 0.0, 4.0, 4.5, 0.5, 3.0, 0.0};
    const double pixel = 3.0;
    const int cols = static_cast<int>(profile.size());
    const Beam beam = makeBeam(cols, 1, pixel);
    const double firstCenter = -0.5 * (cols - 1) * pixel;
    std::vector<double> z = {-1.5, 1.5};

    SlidingWindowRow row;
    row.profile = profile;
    row.firstCenter = firstCenter;
    row.pixel = pixel;
    row.start = firstCenter - 0.5 * pixel;      // one pixel before the first exposed pixel (index 1)
    row.end = firstCenter + (cols - 1.5) * pixel;   // right edge of the last exposed pixel (index cols-2)

    for (int slots : {400, 40}) {
        const double speed = 30.0 / slots * 4.0;   // sweeps the field within a quarter of the slots
        const double mu = SlidingWindow::requiredMuPerSlot(row, speed, slots);
        ASSERT_GT(mu, 0.0);
        std::vector<double> trailing, leading;
        SlidingWindow::sampleTips(row, speed, mu, slots, trailing, leading);

        std::vector<double> delivered(profile.size(), 0.0);
        for (int k = 0; k < slots; ++k) {
            Aperture previous, current, next;
            auto make = [&](int index) {
                const int clamped = std::clamp(index, 0, slots - 1);
                Aperture ap; ap.bankA = {trailing[static_cast<size_t>(clamped)]};
                ap.bankB = {leading[static_cast<size_t>(clamped)]}; return ap; };
            previous = make(k - 1); current = make(k); next = make(k + 1);
            const auto open = ApertureFluence::integratedOpenFractions(beam, &previous, current, &next, z);
            for (size_t i = 0; i < delivered.size(); ++i) delivered[i] += open[i] * mu;
        }
        double error = 0.0, norm = 0.0;
        for (size_t i = 0; i < profile.size(); ++i) {
            error += (delivered[i] - profile[i]) * (delivered[i] - profile[i]);
            norm += profile[i] * profile[i];
        }
        EXPECT_LT(std::sqrt(error / norm), slots == 400 ? 0.01 : 0.05) << slots << " slots";
    }
}

// ============================================================================
// Arc sequencing and timing
// ============================================================================

namespace {

struct ArcFixture {
    std::shared_ptr<PatientData> patient = WaterPhantom::buildCentredPatient(5.0, 80.0, 20.0);
    PlanBuilder::Result built;

    ArcFixture() {
        PlanConfig pc;
        pc.bixelWidth = 10.0;
        pc.deliveryMode = DeliveryMode::VMAT;
        pc.vmat.arcs[0].startAngle = 0;
        pc.vmat.arcs[0].stopAngle = 90;
        pc.vmat.arcs[0].clockwise = true;
        pc.vmat.optimizationSpacingDeg = 30.0;
        pc.vmat.controlPointSpacingDeg = 10.0;
        built = PlanBuilder::build(pc, patient);
    }
};

} // namespace

TEST(VmatPlanTest, PlanBuilderCreatesOneBeamPerSectorAtTheSectorCentre) {
    ArcFixture f;
    ASSERT_NE(f.built.stf, nullptr);
    EXPECT_TRUE(f.built.plan->isVmat());
    ASSERT_EQ(f.built.stf->getCount(), 3u);
    EXPECT_NEAR(f.built.stf->getBeam(0)->getGantryAngle(), 15.0, 1e-9);
    EXPECT_NEAR(f.built.stf->getBeam(2)->getGantryAngle(), 75.0, 1e-9);
}

TEST(ArcSequencerTest, LevelSetUniformFluenceGivesOneSharedApertureAndMuIsConserved) {
    ArcFixture f;
    const Stf& stf = *f.built.stf;
    std::vector<double> weights(stf.getTotalNumOfBixels(), 1.0);

    ArcSequencerOptions levels;
    levels.method = ArcSequencerOptions::Method::LevelSet;
    const auto seq = ArcSequencer::sequence(weights, stf, f.built.plan->getVmatSettings(),
                                            f.built.plan->getMachine(), levels);
    ASSERT_EQ(seq.controlPoints.size(), 9u);

    // One segment per sector (uniform fluence), MU 1 per sector split over its 3 control points.
    EXPECT_NEAR(seq.totalMU, 3.0, 1e-9);
    for (size_t s = 0; s < 3; ++s) {
        const auto& a = seq.controlPoints[3 * s].aperture;
        for (size_t k = 1; k < 3; ++k) {
            const auto& b = seq.controlPoints[3 * s + k].aperture;
            EXPECT_EQ(a.bankA, b.bankA);
            EXPECT_EQ(a.bankB, b.bankB);
        }
        for (size_t k = 0; k < 3; ++k) EXPECT_NEAR(seq.controlPoints[3 * s + k].mu, 1.0 / 3.0, 1e-9);
    }

    // No leaf movement inside a sector, MU is tiny: the gantry runs at full speed.
    EXPECT_EQ(seq.slowedControlPoints, 0);
    EXPECT_NEAR(seq.minGantrySpeed, f.built.plan->getMachine().getConstraints().gantryRotationSpeed[1], 1e-9);
    EXPECT_NEAR(seq.totalTimeSec, 90.0 / 6.0, 1e-6);
}

TEST(ArcSequencerTest, SweepDeliversTheSectorFluenceWithMovingLeaves) {
    ArcFixture f;
    const Stf& stf = *f.built.stf;
    std::vector<double> weights(stf.getTotalNumOfBixels());
    size_t at = 0;
    for (size_t s = 0; s < stf.getCount(); ++s) {
        const Beam& beam = *stf.getBeam(s);
        for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
            const auto& pos = beam.getRay(r)->getRayPosBev();
            weights[at++] = 1.0 + 0.5 * std::sin(0.15 * pos[0]) + 0.3 * std::cos(0.2 * pos[2]);
        }
    }

    const auto& plan = *f.built.plan;
    const auto seq = ArcSequencer::sequence(weights, stf, plan.getVmatSettings(), plan.getMachine());
    ASSERT_EQ(seq.controlPoints.size(), 9u);
    for (const auto& cp : seq.controlPoints) {
        EXPECT_GT(cp.mu, 0.0);
        ASSERT_EQ(cp.aperture.bankA.size(), cp.leafPairBoundariesZ.size() - 1);
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) EXPECT_LE(cp.aperture.bankA[p], cp.aperture.bankB[p]);
    }

    // Delivered fluence (leaves integrated over each slice) against the optimized fluence of each sector.
    size_t offset = 0;
    for (size_t s = 0; s < stf.getCount(); ++s) {
        const Beam& beam = *stf.getBeam(s);
        std::vector<double> delivered(beam.getNumOfRays(), 0.0);
        for (size_t j = 3 * s; j < 3 * s + 3; ++j) {
            const auto& cp = seq.controlPoints[j];
            const Aperture* before = j > 0 ? &seq.controlPoints[j - 1].aperture : nullptr;
            const Aperture* after = j + 1 < seq.controlPoints.size() ? &seq.controlPoints[j + 1].aperture : nullptr;
            const auto open = ApertureFluence::integratedOpenFractions(beam, before, cp.aperture, after,
                                                                       cp.leafPairBoundariesZ);
            for (size_t r = 0; r < open.size(); ++r) delivered[r] += open[r] * cp.mu;
        }
        double error = 0.0, norm = 0.0, got = 0.0, want = 0.0;
        for (size_t r = 0; r < delivered.size(); ++r) {
            const double t = weights[offset + r];
            error += (delivered[r] - t) * (delivered[r] - t);
            norm += t * t;
            got += delivered[r];
            want += t;
        }
        offset += beam.getNumOfRays();
        EXPECT_NEAR(got / want, 1.0, 0.1) << "sector " << s;
        EXPECT_LT(std::sqrt(error / norm), 0.35) << "sector " << s;
    }
}

TEST(ArcSequencerTest, SweepNeedsAnInterdigitatingMlcAndSeveralControlPointsPerSector) {
    ArcFixture f;
    std::vector<double> weights(f.built.stf->getTotalNumOfBixels(), 1.0);
    ArcSequencerOptions sweep;
    sweep.method = ArcSequencerOptions::Method::Sweep;
    auto settings = f.built.plan->getVmatSettings();
    settings.controlPointSpacingDeg = 30.0;   // one control point per sector
    EXPECT_THROW(ArcSequencer::sequence(weights, *f.built.stf, settings, f.built.plan->getMachine(), sweep),
                 std::invalid_argument);
    // Auto falls back to the level decomposition.
    EXPECT_NO_THROW(ArcSequencer::sequence(weights, *f.built.stf, settings, f.built.plan->getMachine()));
}

TEST(SlidingWindowTest, TotalRiseCountsOnlyTheUpwardSteps) {
    EXPECT_DOUBLE_EQ(SlidingWindow::totalRise({0.0, 2.0, 1.0, 4.0, 0.0}), 2.0 + 3.0);
    EXPECT_DOUBLE_EQ(SlidingWindow::totalRise({3.0, 3.0}), 3.0);
    EXPECT_DOUBLE_EQ(SlidingWindow::totalRise({0.0, 0.0}), 0.0);
}

TEST(FluenceMapTest, ToBeamWeightsInvertsFromBeamWeights) {
    const Beam beam = makeBeam(4, 3, 5.0);
    std::vector<double> weights(beam.getNumOfRays());
    for (size_t i = 0; i < weights.size(); ++i) weights[i] = 1.0 + static_cast<double>(i);
    FluenceMap map = FluenceMap::fromBeamWeights(beam, weights, 0);
    map.setValue(1, 2, 100.0);
    std::vector<double> back(weights.size(), -1.0);
    map.toBeamWeights(beam, back, 0);
    for (size_t i = 0; i < weights.size(); ++i)
        EXPECT_DOUBLE_EQ(back[i], i == 1 * 4 + 2 ? 100.0 : weights[i]);
}

TEST(ArcSequencerTest, RiseLimitReducesTheMuOfAModulatedSectorAndKeepsItsMass) {
    ArcFixture f;
    const Stf& stf = *f.built.stf;
    std::vector<double> weights(stf.getTotalNumOfBixels());
    size_t at = 0;
    for (size_t s = 0; s < stf.getCount(); ++s) {
        const Beam& beam = *stf.getBeam(s);
        for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
            const auto& pos = beam.getRay(r)->getRayPosBev();
            const int column = static_cast<int>(std::lround(pos[0] / 10.0));
            weights[at++] = (column % 2 == 0) ? 4.0 : 0.5;   // comb: large upward steps in every row
        }
    }
    const auto& plan = *f.built.plan;
    ArcSequencerOptions raw;
    raw.maxRiseRatio = 0.0;
    raw.fitIterations = 0;
    ArcSequencerOptions limited = raw;
    limited.maxRiseRatio = 1.0;
    const auto a = ArcSequencer::sequence(weights, stf, plan.getVmatSettings(), plan.getMachine(), raw);
    const auto b = ArcSequencer::sequence(weights, stf, plan.getVmatSettings(), plan.getMachine(), limited);
    EXPECT_GT(a.totalMU, 0.0);
    EXPECT_LT(b.totalMU, 0.8 * a.totalMU);
}

TEST(ArcSequencerTest, SweepSlicesAreFastEnoughForTheLeavesAndDoNotExceedTheDoseRate) {
    ArcFixture f;
    std::vector<double> weights(f.built.stf->getTotalNumOfBixels(), 1.0);
    for (size_t i = 0; i < weights.size(); ++i) weights[i] = 1.0 + 0.8 * static_cast<double>((i * 5) % 7) / 7.0;
    const auto& plan = *f.built.plan;
    const auto seq = ArcSequencer::sequence(weights, *f.built.stf, plan.getVmatSettings(), plan.getMachine());
    const auto limits = plan.getMachine().getConstraints();
    for (const auto& cp : seq.controlPoints) {
        EXPECT_LE(cp.leafSpeed, limits.leafSpeed[1] + 1e-6);
        EXPECT_LE(cp.doseRate, limits.monitorUnitRate[1] + 1e-6);
        EXPECT_LE(cp.gantrySpeed, limits.gantryRotationSpeed[1] + 1e-6);
    }
}

TEST(ArcSequencerTest, ConformalAperturesStayOpenAndMoveLittle) {
    ArcFixture f;
    const Stf& stf = *f.built.stf;
    std::vector<double> weights(stf.getTotalNumOfBixels());
    size_t at = 0;
    for (size_t s = 0; s < stf.getCount(); ++s) {
        const Beam& beam = *stf.getBeam(s);
        for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
            const auto& pos = beam.getRay(r)->getRayPosBev();
            weights[at++] = 1.0 + 0.5 * std::sin(0.15 * pos[0]);
        }
    }
    ArcSequencerOptions conformal;
    conformal.method = ArcSequencerOptions::Method::Conformal;
    ArcSequencerOptions sweep;
    sweep.method = ArcSequencerOptions::Method::Sweep;
    const auto& plan = *f.built.plan;
    const auto c = ArcSequencer::sequence(weights, stf, plan.getVmatSettings(), plan.getMachine(), conformal);
    const auto w = ArcSequencer::sequence(weights, stf, plan.getVmatSettings(), plan.getMachine(), sweep);
    ASSERT_EQ(c.controlPoints.size(), 9u);
    for (const auto& cp : c.controlPoints) {
        EXPECT_GT(cp.mu, 0.0);
        bool open = false;
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) open = open || cp.aperture.bankB[p] > cp.aperture.bankA[p];
        EXPECT_TRUE(open);
    }
    EXPECT_LT(c.totalMU, w.totalMU);
    EXPECT_LE(c.totalTimeSec, w.totalTimeSec + 1e-9);
}

TEST(ArcSequencerTest, ResampleInterpolatesTipsAndMuPerDegreeAlongTheArc) {
    ArcSequenceResult source;
    for (int k = 0; k < 2; ++k) {
        DeliveredControlPoint cp;
        cp.geometry.offset = 20.0 + 40.0 * k;
        cp.geometry.width = 40.0;
        cp.leafPairBoundariesZ = {-5.0, 5.0, 15.0};
        cp.aperture.bankA = {0.0 + 10.0 * k, 10.0 * k};
        cp.aperture.bankB = {10.0 + 20.0 * k, 30.0 + 20.0 * k};
        cp.mu = 4.0 + 4.0 * k;   // 0.1 and 0.2 MU per degree
        source.controlPoints.push_back(cp);
    }
    std::vector<ArcControlPoint> geometry(4);
    for (int j = 0; j < 4; ++j) {
        geometry[j].offset = 10.0 + 20.0 * j;
        geometry[j].width = 20.0;
    }
    const auto out = ArcSequencer::resample(source, geometry, MachineConstraints{});
    ASSERT_EQ(out.controlPoints.size(), 4u);
    EXPECT_NEAR(out.controlPoints[0].mu, 0.1 * 20.0, 1e-12);          // before the first source point: held
    EXPECT_NEAR(out.controlPoints[1].mu, 0.125 * 20.0, 1e-12);        // a quarter of the way
    EXPECT_NEAR(out.controlPoints[2].mu, 0.175 * 20.0, 1e-12);        // three quarters
    EXPECT_NEAR(out.controlPoints[3].mu, 0.2 * 20.0, 1e-12);          // after the last: held
    EXPECT_NEAR(out.controlPoints[1].aperture.bankB[0], 10.0 + 0.25 * 20.0, 1e-12);
    EXPECT_NEAR(out.controlPoints[2].aperture.bankA[1], 7.5, 1e-12);
    EXPECT_EQ(out.controlPoints[0].leafPairBoundariesZ, source.controlPoints[0].leafPairBoundariesZ);
}

TEST(ArcSequencerTest, RejectsAnStfThatDoesNotMatchTheArcs) {
    ArcFixture f;
    VmatSettings other = f.built.plan->getVmatSettings();
    other.optimizationSpacingDeg = 10.0;   // 9 sectors instead of 3
    std::vector<double> weights(f.built.stf->getTotalNumOfBixels(), 1.0);
    EXPECT_THROW(ArcSequencer::sequence(weights, *f.built.stf, other, f.built.plan->getMachine()),
                 std::invalid_argument);
}

TEST(ArcSequencerTest, FastLeafMotionSlowsTheGantryAndHighMuSlowsItToo) {
    ArcFixture f;
    std::vector<double> weights(f.built.stf->getTotalNumOfBixels(), 1.0);
    ArcSequencerOptions levels;
    levels.method = ArcSequencerOptions::Method::LevelSet;   // static apertures: no sweep travel inside a sector
    auto seq = ArcSequencer::sequence(weights, *f.built.stf, f.built.plan->getVmatSettings(),
                                      f.built.plan->getMachine(), levels);
    const MachineConstraints limits = f.built.plan->getMachine().getConstraints();

    // Throw the leaves 200 mm in one slice: needs 200/60 s instead of 10/6 s.
    for (auto& v : seq.controlPoints[4].aperture.bankA) v += 200.0;
    for (auto& v : seq.controlPoints[4].aperture.bankB) v += 200.0;
    ArcSequencer::computeTiming(seq, limits);
    EXPECT_GE(seq.slowedControlPoints, 1);
    EXPECT_NEAR(seq.controlPoints[4].timeSec, 200.0 / limits.leafSpeed[1], 1e-9);
    EXPECT_LT(seq.controlPoints[4].gantrySpeed, limits.gantryRotationSpeed[1]);
    EXPECT_NEAR(seq.maxLeafSpeed, limits.leafSpeed[1], 1e-9);

    // 100 MU in one slice at 10 MU/s needs 10 s.
    seq.controlPoints[0].mu = 100.0;
    ArcSequencer::computeTiming(seq, limits);
    EXPECT_NEAR(seq.controlPoints[0].timeSec, 100.0 / limits.monitorUnitRate[1], 1e-9);
    EXPECT_NEAR(seq.controlPoints[0].doseRate, limits.monitorUnitRate[1], 1e-9);
}

TEST(ArcSequencerTest, LeafSpeedLimitKeepsTheGantryAtFullSpeed) {
    ArcFixture f;
    std::vector<double> weights(f.built.stf->getTotalNumOfBixels());
    for (size_t i = 0; i < weights.size(); ++i) weights[i] = 0.01 * (1.0 + static_cast<double>((i * 7) % 11));

    ArcSequencerOptions opt;
    opt.enforceLeafSpeed = true;
    const auto seq = ArcSequencer::sequence(weights, *f.built.stf, f.built.plan->getVmatSettings(),
                                            f.built.plan->getMachine(), opt);
    const MachineConstraints limits = f.built.plan->getMachine().getConstraints();

    EXPECT_LE(seq.maxLeafSpeed, limits.leafSpeed[1]);
    EXPECT_EQ(seq.slowedControlPoints, 0);
    EXPECT_NEAR(seq.minGantrySpeed, limits.gantryRotationSpeed[1], 1e-9);
}

TEST(ArcSequencerTest, ClosedLeavesMovePhysicallyAndArcsDoNotShareTransitions) {
    ArcSequenceResult sequence;
    DeliveredControlPoint first;
    first.geometry.width = 4.0;
    first.leafPairBoundariesZ = {-5.0, 5.0};
    first.aperture.bankA = {0.0};
    first.aperture.bankB = {0.0};
    first.mu = 10.0;
    auto second = first;
    second.aperture.bankA = {120.0};
    second.aperture.bankB = {140.0};
    sequence.controlPoints = {first, second};
    MachineConstraints limits;
    ArcSequencer::computeTiming(sequence, limits);
    EXPECT_NEAR(sequence.controlPoints[1].timeSec, 140.0 / limits.leafSpeed[1], 1e-9);
    EXPECT_NEAR(sequence.controlPoints[1].leafSpeed, limits.leafSpeed[1], 1e-9);
    sequence.controlPoints[1].geometry.arcIndex = 1;
    ArcSequencer::computeTiming(sequence, limits);
    EXPECT_DOUBLE_EQ(sequence.controlPoints[1].leafSpeed, 0.0);
}

// ============================================================================
// Full arc pipeline on a water cube
// ============================================================================

namespace {

BuiltObjectives ptvObjective(const PatientData& patient, const Grid& doseGrid, double dose) {
    ObjectiveSpec spec;
    spec.structurePattern = "PTV";
    spec.typeIdx = 3;           // MinDVH
    spec.doseValue = dose;
    spec.weight = 100.0;
    spec.volumePct = 98.0;
    return ObjectiveBuilder::buildFromSpecs({spec}, patient, patient.getGrid(), doseGrid);
}

double ptvMeanDose(const ArcPipelineResult& r) {
    for (const auto& s : r.stats)
        if (s.name == "PTV") return s.meanDose;
    return -1.0;
}

} // namespace

TEST(ArcPipelineTest, WeightOptimizationBringsTheTargetToThePrescription) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});

    ArcOptimizationConfig cfg;
    cfg.optimizeLeaves = false;
    cfg.computeBackend = DoseComputeBackend::CPU;
    cfg.base.targetDose = 60.0;
    cfg.weightIterations = 100;
    cfg.base.ntoEnabled = false;
    cfg.muSmoothingWeight = 0.0;

    // Calibrate: dose per unit fluence weight with the sequenced MU delivered as is.
    std::vector<double> unit(f.built.stf->getTotalNumOfBixels(), 1.0);
    cfg.optimizeWeights = false;
    const auto calib = ArcPipeline::run(unit, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                        BuiltObjectives{}, cfg);
    const double perUnit = ptvMeanDose(calib);
    ASSERT_GT(perUnit, 0.0);

    // Start at 60 % of the prescription; the optimizer must raise the MU.
    std::vector<double> weights(unit.size(), 0.6 * 60.0 / perUnit);
    cfg.optimizeWeights = true;
    const auto r = ArcPipeline::run(weights, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                    ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    EXPECT_GT(r.iterations, 0);
    EXPECT_GT(r.sequence.totalMU, calib.sequence.totalMU * 0.6 * 60.0 / perUnit);  // MU went up
    EXPECT_NEAR(ptvMeanDose(r), 60.0, 0.1 * 60.0);

    ASSERT_NE(r.dose, nullptr);
    EXPECT_GT(r.dose->getMax(), 0.0);
    EXPECT_EQ(r.controlPointStf->getCount(), 9u);
    EXPECT_TRUE(std::isfinite(r.sequence.totalTimeSec));
    EXPECT_GE(r.sequence.totalTimeSec, 90.0 / 6.0 - 1e-9);
}

TEST(ArcPipelineTest, LeafRefinementDecreasesDoseObjectiveAndRespectsRates) {
    ArcFixture fixture;
    const Grid grid = Grid::createDoseGrid(fixture.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig config;
    config.optimizeWeights = false;
    config.base.ntoEnabled = false;
    config.base.targetDose = 60.0;
    config.leafIterations = 5;
    std::vector<double> weights(fixture.built.stf->getTotalNumOfBixels(), 30.0);
    ObjectiveSpec spec;
    spec.structurePattern = "PTV";
    spec.typeIdx = 0;
    spec.doseValue = 60.0;
    spec.weight = 100.0;
    auto objectives = ObjectiveBuilder::buildFromSpecs({spec}, *fixture.patient, fixture.patient->getGrid(), grid);
    const auto result = ArcPipeline::run(weights, *fixture.built.stf, *fixture.built.plan,
        *fixture.patient, grid, std::move(objectives), config);
    EXPECT_GT(result.leafIterations, 0);
    EXPECT_LT(result.finalLeafObjective, result.initialLeafObjective);
    const auto& limits = fixture.built.plan->getMachine().getConstraints();
    for (const auto& cp : result.sequence.controlPoints) {
        EXPECT_LE(cp.leafSpeed, limits.leafSpeed[1] + 1e-8);
        EXPECT_LE(cp.gantrySpeed, limits.gantryRotationSpeed[1] + 1e-8);
        EXPECT_LE(cp.doseRate, limits.monitorUnitRate[1] + 1e-8);
        if (cp.mu > 0.0) EXPECT_GE(cp.doseRate, limits.monitorUnitRate[0] - 1e-8);
        for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
            EXPECT_LE(cp.aperture.bankA[pair], cp.aperture.bankB[pair]);
            const double resolution = fixture.built.plan->getMachine().getGeometry().leafPositionResolution;
            if (resolution > 0.0) {
                EXPECT_NEAR(cp.aperture.bankA[pair] / resolution,
                    std::round(cp.aperture.bankA[pair] / resolution), 1e-9);
                EXPECT_NEAR(cp.aperture.bankB[pair] / resolution,
                    std::round(cp.aperture.bankB[pair] / resolution), 1e-9);
            }
        }
    }
    for (size_t index = 1; index < result.sequence.controlPoints.size(); ++index) {
        const auto& cp = result.sequence.controlPoints[index];
        const auto& previous = result.sequence.controlPoints[index - 1];
        if (cp.geometry.arcIndex != previous.geometry.arcIndex) continue;
        ASSERT_EQ(cp.leafPairBoundariesZ, previous.leafPairBoundariesZ);
        for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
            EXPECT_LE(std::abs(cp.aperture.bankA[pair] - previous.aperture.bankA[pair]),
                limits.leafSpeed[1] * cp.timeSec + 1e-8);
            EXPECT_LE(std::abs(cp.aperture.bankB[pair] - previous.aperture.bankB[pair]),
                limits.leafSpeed[1] * cp.timeSec + 1e-8);
        }
    }
}


TEST(ArcPipelineTest, MuPenaltyLowersTheTotalMuWhileTheTargetStaysCovered) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig cfg;
    cfg.computeBackend = DoseComputeBackend::CPU;
    cfg.base.targetDose = 60.0;
    cfg.base.ntoEnabled = false;
    cfg.muSmoothingWeight = 0.0;
    cfg.leafIterations = 30;
    cfg.maxWeightFactor = 3.0;

    std::vector<double> unit(f.built.stf->getTotalNumOfBixels(), 1.0);
    ArcOptimizationConfig calibration = cfg;
    calibration.optimizeWeights = false;
    calibration.optimizeLeaves = false;
    const double perUnit = ptvMeanDose(ArcPipeline::run(unit, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                                         BuiltObjectives{}, calibration));
    ASSERT_GT(perUnit, 0.0);
    std::vector<double> weights(unit.size(), 60.0 / perUnit);

    cfg.muPenalty = 0.0;
    const auto free = ArcPipeline::run(weights, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                       ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    cfg.muPenalty = 1.0;
    const auto penalized = ArcPipeline::run(weights, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                            ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    EXPECT_LT(penalized.sequence.totalMU, free.sequence.totalMU);
    EXPECT_GT(ptvMeanDose(penalized), 0.8 * 60.0);
}


TEST(ArcPipelineTest, ProgressiveDaoRefinesThroughSeveralResolutionsAndRespectsTheLimits) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig cfg;
    cfg.computeBackend = DoseComputeBackend::CPU;
    cfg.base.targetDose = 60.0;
    cfg.base.ntoEnabled = false;
    cfg.strategy = ArcStrategy::ProgressiveDao;
    cfg.daoCoarsestSpacingDeg = 40.0;   // 40, 20, 10 deg control points: 2, 5, 9
    cfg.daoIterations = 15;

    std::vector<double> unit(f.built.stf->getTotalNumOfBixels(), 1.0);
    ArcOptimizationConfig calibration = cfg;
    calibration.strategy = ArcStrategy::SweepRefine;
    calibration.optimizeWeights = false;
    calibration.optimizeLeaves = false;
    const double perUnit = ptvMeanDose(ArcPipeline::run(unit, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                                         BuiltObjectives{}, calibration));
    ASSERT_GT(perUnit, 0.0);
    std::vector<double> weights(unit.size(), 60.0 / perUnit);

    const auto r = ArcPipeline::run(weights, *f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                    ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    ASSERT_EQ(r.sequence.controlPoints.size(), 9u);
    EXPECT_EQ(r.controlPointStf->getCount(), 9u);
    EXPECT_GT(r.leafIterations, 0);
    EXPECT_LT(r.finalLeafObjective, r.initialLeafObjective);
    EXPECT_NEAR(ptvMeanDose(r), 60.0, 0.25 * 60.0);

    const auto limits = f.built.plan->getMachine().getConstraints();
    for (const auto& cp : r.sequence.controlPoints) {
        EXPECT_LE(cp.leafSpeed, limits.leafSpeed[1] + 1e-6);
        EXPECT_LE(cp.doseRate, limits.monitorUnitRate[1] + 1e-6);
        EXPECT_LE(cp.gantrySpeed, limits.gantryRotationSpeed[1] + 1e-6);
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) EXPECT_LE(cp.aperture.bankA[p], cp.aperture.bankB[p]);
    }
}

} // namespace optirad::tests
