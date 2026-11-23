#include <gtest/gtest.h>

#include "core/ArcPlan.hpp"
#include "core/Beam.hpp"
#include "core/FluenceMap.hpp"
#include "core/Ray.hpp"
#include "core/workflow/ArcPipeline.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "optimization/OptimizationLogging.hpp"
#include "optimization/ObjectiveProtocol.hpp"
#include "dose/engines/PencilBeamEngine.hpp"
#include "sequencing/ApertureFluence.hpp"
#include "sequencing/ArcSequencer.hpp"
#include "validation/WaterPhantom.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>

namespace optirad::tests {

// ============================================================================
// Arc geometry
// ============================================================================

TEST(ArcGeometryTest, PlanningDefaultsCoverAFullArcWithRequestedResolution) {
    const PlanConfig config;
    EXPECT_EQ(config.machineName, "HD120");
    EXPECT_EQ(config.deliveryMode, DeliveryMode::VMAT);
    EXPECT_DOUBLE_EQ(config.bixelWidth, 3.0);
    const auto& settings = config.vmat;
    ASSERT_EQ(settings.arcs.size(), 1u);
    EXPECT_DOUBLE_EQ(settings.arcs[0].startAngle, 0.0);
    EXPECT_DOUBLE_EQ(settings.arcs[0].stopAngle, 360.0);
    EXPECT_TRUE(settings.arcs[0].clockwise);
    EXPECT_DOUBLE_EQ(settings.arcs[0].span(), 360.0);
    EXPECT_DOUBLE_EQ(settings.optimizationSpacingDeg, 10.0);
    EXPECT_DOUBLE_EQ(settings.controlPointSpacingDeg, 2.0);

    const auto sectors = ArcGeometry::sectors(settings);
    const auto cps = ArcGeometry::controlPoints(settings);
    ASSERT_EQ(sectors.size(), 36u);
    ASSERT_EQ(cps.size(), 180u);
    double totalWidth = 0.0;
    for (size_t i = 0; i < cps.size(); ++i) {
        EXPECT_DOUBLE_EQ(cps[i].width, 2.0);
        EXPECT_DOUBLE_EQ(cps[i].gantryAngle, 1.0 + 2.0 * i);
        EXPECT_EQ(cps[i].sectorIndex, i / 5);
        totalWidth += cps[i].width;
    }
    EXPECT_DOUBLE_EQ(totalWidth, 360.0);
}

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

// ============================================================================
// Arc sequencing and timing
// ============================================================================

namespace {

struct ArcFixture {
    std::shared_ptr<PatientData> patient = WaterPhantom::buildCentredPatient(5.0, 80.0, 20.0);
    PlanBuilder::Result built;

    ArcFixture() {
        PlanConfig pc;
        pc.machineName = "Generic";
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

TEST(ArcSequencerTest, ConformalSeedOpensEveryControlPointWithUnitMu) {
    ArcFixture f;
    const auto& plan = *f.built.plan;
    const auto seed = ArcSequencer::conformalSeed(*f.built.stf, plan.getVmatSettings(), plan.getMachine());
    ASSERT_EQ(seed.controlPoints.size(), 9u);
    for (const auto& cp : seed.controlPoints) {
        EXPECT_GT(cp.mu, 0.0);
        bool open = false;
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) {
            EXPECT_LE(cp.aperture.bankA[p], cp.aperture.bankB[p]);
            open = open || cp.aperture.bankB[p] > cp.aperture.bankA[p];
        }
        EXPECT_TRUE(open);
    }
    EXPECT_TRUE(std::isfinite(seed.totalTimeSec));
    EXPECT_GE(seed.totalTimeSec, 90.0 / 6.0 - 1e-9);
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
    EXPECT_THROW(ArcSequencer::conformalSeed(*f.built.stf, other, f.built.plan->getMachine()),
                 std::invalid_argument);
}

TEST(ArcSequencerTest, FastLeafMotionSlowsTheGantryAndHighMuSlowsItToo) {
    ArcFixture f;
    auto seq = ArcSequencer::conformalSeed(*f.built.stf, f.built.plan->getVmatSettings(),
                                           f.built.plan->getMachine());
    for (auto& cp : seq.controlPoints) cp.aperture = seq.controlPoints[4].aperture;   // static leaves
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

size_t countOccurrences(const std::string& text, const std::string& needle) {
    size_t count = 0, position = 0;
    while ((position = text.find(needle, position)) != std::string::npos) {
        ++count;
        position += needle.size();
    }
    return count;
}

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

TEST(ArcPipelineTest, DaoDecreasesTheObjectiveAndRespectsTheMachineLimits) {
    ArcFixture fixture;
    const Grid grid = Grid::createDoseGrid(fixture.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig config;
    config.computeBackend = DoseComputeBackend::CPU;
    config.base.nto.enabled = false;
    config.base.targetDose = 60.0;
    config.coarsestSpacingDeg = 10.0;
    config.iterationsPerResolution = 5;
    ObjectiveSpec spec;
    spec.structurePattern = "PTV";
    spec.typeIdx = 0;
    spec.doseValue = 60.0;
    spec.weight = 100.0;
    auto objectives = ObjectiveBuilder::buildFromSpecs({spec}, *fixture.patient, fixture.patient->getGrid(), grid);
    testing::internal::CaptureStdout();
    const auto result = ArcPipeline::run(*fixture.built.stf, *fixture.built.plan,
        *fixture.patient, grid, std::move(objectives), config);
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_EQ(countOccurrences(log, "weighted cost breakdown total="), 2u);
    EXPECT_NE(log.find("DAO initial: weighted cost breakdown"), std::string::npos);
    EXPECT_NE(log.find("DAO final (after feasibility repair): weighted cost breakdown"), std::string::npos);
    EXPECT_GT(result.leafIterations, 0);
    EXPECT_LT(result.finalLeafObjective, result.initialLeafObjective);
    // `mu` are course bixel weights; the dose-rate limits are converted to them (1 = 100 MU x (SAD/1 m)^2 / fx).
    const auto& plan = *fixture.built.plan;
    EXPECT_NEAR(ArcPipeline::machineMuPerWeight(plan), 100.0 * std::pow(plan.getMachine().getSAD() / 1000.0, 2) /
                plan.getNumOfFractions(), 1e-12);
    const auto limits = ArcPipeline::deliveryConstraints(plan);
    EXPECT_NEAR(limits.monitorUnitRate[1] * ArcPipeline::machineMuPerWeight(plan),
                plan.getMachine().getConstraints().monitorUnitRate[1], 1e-9);
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
    cfg.base.nto.enabled = false;
    cfg.muSmoothing = 0.0;
    cfg.coarsestSpacingDeg = 20.0;
    cfg.iterationsPerResolution = 30;

    cfg.muPenalty = 0.0;
    const auto free = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                       ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    cfg.muPenalty = 1.0;
    const auto penalized = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
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
    cfg.base.nto.enabled = false;
    cfg.coarsestSpacingDeg = 40.0;   // 40, 20, 10 deg control points: 2, 5, 9
    cfg.iterationsPerResolution = 15;

    const auto r = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                    ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    ASSERT_EQ(r.sequence.controlPoints.size(), 9u);
    EXPECT_EQ(r.controlPointStf->getCount(), 9u);
    EXPECT_GT(r.leafIterations, 0);
    EXPECT_LT(r.finalLeafObjective, r.initialLeafObjective);
    EXPECT_NEAR(ptvMeanDose(r), 60.0, 0.25 * 60.0);

    const auto limits = ArcPipeline::deliveryConstraints(*f.built.plan);
    for (const auto& cp : r.sequence.controlPoints) {
        EXPECT_LE(cp.leafSpeed, limits.leafSpeed[1] + 1e-6);
        EXPECT_LE(cp.doseRate, limits.monitorUnitRate[1] + 1e-6);
        EXPECT_LE(cp.gantrySpeed, limits.gantryRotationSpeed[1] + 1e-6);
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) EXPECT_LE(cp.aperture.bankA[p], cp.aperture.bankB[p]);
    }
}

TEST(ArcPipelineTest, ProgressiveDaoReportsProgressAndCalibratesItsOwnMu) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig cfg;
    cfg.computeBackend = DoseComputeBackend::CPU;
    cfg.base.targetDose = 60.0;
    cfg.base.nto.enabled = false;
    cfg.coarsestSpacingDeg = 40.0;
    cfg.iterationsPerResolution = 15;

    auto objectives = ptvObjective(*f.patient, doseGrid, 60.0);
    std::vector<ArcProgressSnapshot> snapshots;
    testing::internal::CaptureStdout();
    const auto r = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                    std::move(objectives), cfg, nullptr,
                                    [&](const ArcProgressSnapshot& sample) { snapshots.push_back(sample); });
    const auto log = testing::internal::GetCapturedStdout();
    ASSERT_EQ(snapshots.size(), r.progressHistory.size());
    ASSERT_GE(snapshots.size(), 6u);
    int completedStages = 0;
    int lastStage = 0, lastIteration = 0;
    for (size_t sampleIndex = 0; sampleIndex < snapshots.size(); ++sampleIndex) {
        const auto& sample = snapshots[sampleIndex];
        EXPECT_EQ(sample.stages, 3);
        EXPECT_EQ(sample.stage, r.progressHistory[sampleIndex].stage);
        if (sample.stage != lastStage) {
            EXPECT_EQ(sample.stage, lastStage + 1);
            EXPECT_EQ(sample.progress.iteration, 0);
            lastStage = sample.stage;
            lastIteration = 0;
        }
        EXPECT_GE(sample.progress.iteration, lastIteration);
        lastIteration = sample.progress.iteration;
        double sum = 0.0;
        for (const auto& term : sample.breakdown.contributions) sum += term.second;
        EXPECT_DOUBLE_EQ(sum, sample.breakdown.total);
        EXPECT_DOUBLE_EQ(sample.progress.objective, sum);
        EXPECT_GT(sample.progress.totalMu, 0.0);
        if (sample.stageFinished) ++completedStages;
        else if (sample.progress.acceptedSteps > 0)
            EXPECT_TRUE(shouldLogOptimizationProgress(sample.progress.acceptedSteps, 10));
    }
    EXPECT_EQ(completedStages, 3);
    EXPECT_TRUE(snapshots.back().stageFinished);
    EXPECT_NEAR(snapshots.back().breakdown.total, r.finalLeafObjective,
                1e-8 * std::max(1.0, std::abs(r.finalLeafObjective)));
    testing::internal::CaptureStdout();
    const auto baseline = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                           ptvObjective(*f.patient, doseGrid, 60.0), cfg);
    testing::internal::GetCapturedStdout();
    EXPECT_TRUE(baseline.progressHistory.empty());
    EXPECT_EQ(baseline.leafIterations, r.leafIterations);
    EXPECT_NEAR(baseline.finalLeafObjective, r.finalLeafObjective,
                1e-8 * std::max(1.0, std::abs(r.finalLeafObjective)));
    EXPECT_NEAR(baseline.sequence.totalMU, r.sequence.totalMU, 1e-8);
    for (const auto* expected : {
            "starting progressive DAO", "weighted dose objectives (soft penalties)",
            "DAO stage 1/3: creating conformal seed", "DAO stage 2/3: interpolating",
            "DAO stage 3/3: interpolating", "DAO seed: calibrating uniform MU",
            "DAO leaf refinement: dose/objective backend=CPU", "DAO settings: iterations=15",
            "Bounds: tips=", "MU total=", "DAO initial: weighted cost breakdown",
            "=== DAO Optimizer (Memory: 10, Max Iter: 15) ===", "Objective", "GradInf", "Time(s)",
            "MU smoothing: cost=",
            "Total MU: cost=", "Leaf order: cost=", "Leaf speed: cost=",
            "DAO final (after feasibility repair)", "DAO stopped:", "Restarts ", "(dir/LS)",
            "DoseCost ", "computing delivered dose",
            "ArcPipeline: complete:"})
        EXPECT_NE(log.find(expected), std::string::npos) << expected;
    EXPECT_EQ(countOccurrences(log, "weighted cost breakdown total="), 2u);
    EXPECT_EQ(countOccurrences(log, "DAO stopped:"), 3u);
    EXPECT_EQ(log.find("DAO accepted step"), std::string::npos);
    EXPECT_EQ(log.find("pre-step gradient inf-norm="), std::string::npos);
    EXPECT_EQ(log.find("DAO iteration "), std::string::npos);
    std::istringstream lines(log);
    int progressLines = 0;
    bool compactBlock = false;
    const std::string prefix = "[INFO] ";
    for (std::string line; std::getline(lines, line);) {
        if (line.find(prefix + "DAO settings:") == 0 || line.find(prefix + "DAO stopped:") == 0)
            compactBlock = true;
        else if (line.find("  ") != 0)
            compactBlock = false;
        if (compactBlock) EXPECT_LE(line.size(), 100u) << line;
        if (line.find(prefix) != 0) continue;
        std::istringstream row(line.substr(prefix.size()));
        DaoProgressInfo info;
        if (!(row >> info.iteration >> info.acceptedSteps >> info.objective >> info.preStepGradInf
                  >> info.stepSize >> info.lsIters >> info.totalMu >> info.elapsedSec)) continue;
        EXPECT_EQ(info.acceptedSteps, 1);
        EXPECT_GE(info.iteration, info.acceptedSteps);
        EXPECT_GE(info.objective, 0.0);
        EXPECT_GE(info.preStepGradInf, 0.0);
        EXPECT_GT(info.stepSize, 0.0);
        EXPECT_GE(info.lsIters, 1);
        EXPECT_LE(info.lsIters, 20);
        EXPECT_GT(info.totalMu, 0.0);
        EXPECT_GE(info.elapsedSec, 0.0);
        EXPECT_LE(line.size(), 100u);
        ++progressLines;
    }
    EXPECT_GE(progressLines, 1);
    EXPECT_LE(progressLines, 3);
    EXPECT_EQ(countOccurrences(log, "=== DAO Optimizer"), static_cast<size_t>(progressLines));
    const std::string finalField = "DAO final (after feasibility repair): weighted cost breakdown total=";
    const auto finalPosition = log.rfind(finalField);
    ASSERT_NE(finalPosition, std::string::npos);
    EXPECT_NEAR(std::stod(log.substr(finalPosition + finalField.size())), r.finalLeafObjective,
                2e-6 * std::max(1.0, std::abs(r.finalLeafObjective)));
    ASSERT_EQ(r.sequence.controlPoints.size(), 9u);
    EXPECT_GT(r.sequence.totalMU, 0.0);
    EXPECT_NEAR(ptvMeanDose(r), 60.0, 0.25 * 60.0);
    for (const auto& cp : r.sequence.controlPoints)
        for (size_t p = 0; p < cp.aperture.bankA.size(); ++p) EXPECT_LE(cp.aperture.bankA[p], cp.aperture.bankB[p]);
}

TEST(ArcPipelineTest, APrescriptionIsRequiredAndNeverAssumed) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    EXPECT_EQ(OptimizationConfig{}.targetDose, 0.0);
    const auto protocol = ObjectiveProtocol::lungIMRT(0.0);
    EXPECT_THROW(ObjectiveBuilder::build(protocol, *f.patient, f.patient->getGrid(), doseGrid, 0.0),
                 std::invalid_argument);
    testing::internal::CaptureStdout();
    const auto built = ObjectiveBuilder::build(protocol, *f.patient, f.patient->getGrid(), doseGrid, 60.0);
    testing::internal::GetCapturedStdout();
    ASSERT_EQ(built.ptrs.size(), 3u);   // target template: Lower D98%, Upper 0% (105%), Squared Deviation
    EXPECT_NE(built.ptrs[1]->describe().find("63"), std::string::npos) << built.ptrs[1]->describe();
    ArcOptimizationConfig cfg;          // no targetDose set
    EXPECT_THROW(ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                  ptvObjective(*f.patient, doseGrid, 60.0), cfg),
                 std::invalid_argument);
}

TEST(ArcPipelineTest, DaoAddsTheNormalTissueObjectiveToTheCostBreakdown) {
    ArcFixture f;
    const Grid doseGrid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    ArcOptimizationConfig cfg;
    cfg.computeBackend = DoseComputeBackend::CPU;
    cfg.base.targetDose = 60.0;
    cfg.coarsestSpacingDeg = 10.0;
    cfg.iterationsPerResolution = 5;
    std::vector<ArcProgressSnapshot> snapshots;
    testing::internal::CaptureStdout();
    const auto r = ArcPipeline::run(*f.built.stf, *f.built.plan, *f.patient, doseGrid,
                                    ptvObjective(*f.patient, doseGrid, 60.0), cfg, nullptr,
                                    [&](const ArcProgressSnapshot& sample) { snapshots.push_back(sample); });
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_NE(log.find("ObjectiveBuilder: added NTO ("), std::string::npos) << log;
    ASSERT_FALSE(snapshots.empty());
    bool hasNto = false;
    for (const auto& term : snapshots.back().breakdown.contributions)
        hasNto = hasNto || term.first.find("NTO (") != std::string::npos;
    EXPECT_TRUE(hasNto);
    EXPECT_EQ(log.find("hotspot"), std::string::npos);
    EXPECT_GT(r.sequence.totalMU, 0.0);
}

TEST(ArcPipelineTest, CoverAperturesAddsRaysWhereAnApertureOpensOutsideTheTarget) {
    ArcFixture f;
    const auto& plan = *f.built.plan;
    auto seq = ArcSequencer::conformalSeed(*f.built.stf, plan.getVmatSettings(), plan.getMachine());
    std::vector<ArcControlPoint> geometry;
    for (const auto& cp : seq.controlPoints) geometry.push_back(cp.geometry);
    auto cpStf = ArcPipeline::buildControlPointStf(plan, *f.patient, geometry);
    ArcPipeline::coverApertures(*cpStf, seq);
    const size_t covered = cpStf->getTotalNumOfRays();
    const size_t beam0Rays = cpStf->getBeam(0)->getNumOfRays();

    auto& aperture = seq.controlPoints[4].aperture;
    for (size_t p = 0; p < aperture.bankA.size(); ++p)
        if (aperture.bankB[p] > aperture.bankA[p]) {
            aperture.bankA[p] = -60.0;
            aperture.bankB[p] = 60.0;
        }
    ArcPipeline::coverApertures(*cpStf, seq);
    const size_t widened = cpStf->getTotalNumOfRays();
    EXPECT_GT(widened, covered);
    EXPECT_EQ(cpStf->getBeam(0)->getNumOfRays(), beam0Rays);
    EXPECT_GT(cpStf->getBeam(4)->getNumOfRays(), beam0Rays);
    EXPECT_GT(cpStf->getBeam(3)->getNumOfRays(), beam0Rays);   // leaf-motion neighbour
    ArcPipeline::coverApertures(*cpStf, seq);
    EXPECT_EQ(cpStf->getTotalNumOfRays(), widened);
}

TEST(ArcPipelineTest, SequenceDoseDeliversAGivenSequenceWithoutOptimization) {
    ArcFixture f;
    const auto& plan = *f.built.plan;
    const Grid grid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    const auto seq = ArcSequencer::conformalSeed(*f.built.stf, plan.getVmatSettings(), plan.getMachine());
    ArcOptimizationConfig config;
    config.computeBackend = DoseComputeBackend::CPU;
    config.base.targetDose = 60.0;
    testing::internal::CaptureStdout();
    const auto r = ArcPipeline::computeSequenceDose(seq, plan, *f.patient, grid, config);
    testing::internal::GetCapturedStdout();
    ASSERT_NE(r.dose, nullptr);
    ASSERT_NE(r.controlPointStf, nullptr);
    EXPECT_EQ(r.controlPointStf->getCount(), seq.controlPoints.size());
    EXPECT_EQ(r.sequence.controlPoints.size(), seq.controlPoints.size());
    EXPECT_DOUBLE_EQ(r.sequence.totalMU, seq.totalMU);
    EXPECT_TRUE(r.progressHistory.empty());
    EXPECT_TRUE(r.fromGivenSequence);
    EXPECT_GT(ptvMeanDose(r), 0.0);

    ArcSequenceResult empty;
    EXPECT_THROW(ArcPipeline::computeSequenceDose(empty, plan, *f.patient, grid, config), std::exception);
}

TEST(ArcPipelineTest, DeliveredDoseHonorsTheBackend) {
    ArcFixture f;
    const Grid grid = Grid::createDoseGrid(f.patient->getGrid(), {5.0, 5.0, 5.0});
    {
        for (const auto backend : {DoseComputeBackend::CPU, DoseComputeBackend::Auto, DoseComputeBackend::CUDA}) {
            SCOPED_TRACE(static_cast<int>(backend));
            ArcOptimizationConfig config;
            config.computeBackend = backend;
            config.coarsestSpacingDeg = 10.0;
            config.iterationsPerResolution = 1;
            config.base.nto.enabled = false;
            config.base.targetDose = 60.0;
            // The selected backend must override stale engine flags in both directions.
            config.engineOptions.gpuRayTracing = backend == DoseComputeBackend::CPU;
            config.engineOptions.requireGpuRayTracing = backend == DoseComputeBackend::CPU;
            auto objectives = ptvObjective(*f.patient, grid, 60.0);
#ifndef OPTIRAD_HAS_CUDA
            if (backend == DoseComputeBackend::CUDA) {
                EXPECT_THROW(ArcPipeline::run(*f.built.stf, *f.built.plan,
                    *f.patient, grid, std::move(objectives), config), std::runtime_error);
                continue;
            }
#endif
            ArcPipelineResult result;
            testing::internal::CaptureStdout();
            EXPECT_NO_THROW(result = ArcPipeline::run(*f.built.stf, *f.built.plan,
                *f.patient, grid, std::move(objectives), config));
            const std::string log = testing::internal::GetCapturedStdout();
            ASSERT_NE(result.dose, nullptr) << log;
            const auto deliveredStart = log.find("ArcPipeline: computing delivered dose");
            ASSERT_NE(deliveredStart, std::string::npos) << log;
            const std::string deliveredLog = log.substr(deliveredStart);
            bool expectCuda = false;
#ifdef OPTIRAD_HAS_CUDA
            expectCuda = backend != DoseComputeBackend::CPU;
#endif
            EXPECT_EQ(deliveredLog.find("CUDA direct-dose deposition enabled") != std::string::npos,
                expectCuda) << deliveredLog;
            EXPECT_EQ(deliveredLog.find("CUDA radiological-depth tracing enabled") != std::string::npos,
                expectCuda) << deliveredLog;

            PencilBeamEngine cpuEngine;
            const auto finalWeights = ArcPipeline::bixelWeights(*result.controlPointStf, result.sequence,
                config.integrateLeafMotion);
            const auto reference = cpuEngine.calculateDoseDirect(*f.built.plan, *result.controlPointStf,
                *f.patient, grid, finalWeights);
            const double peak = reference.getMax();
            ASSERT_GT(peak, 0.0);
            for (size_t voxel = 0; voxel < grid.getNumVoxels(); ++voxel)
                EXPECT_NEAR(result.dose->data()[voxel], reference.data()[voxel], 2e-6 * peak);
        }
    }
}

} // namespace optirad::tests
