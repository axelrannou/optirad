#include <gtest/gtest.h>
#include "core/workflow/ImportedPlanConverter.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "core/Machine.hpp"
#include "core/PatientData.hpp"
#include "core/Plan.hpp"

namespace optirad::tests {

namespace {

Plan planWith(double sad, int fractions) {
    Plan plan;
    Machine machine;
    MachineMeta meta;
    meta.name = "Test";
    meta.SAD = sad;
    machine.setMeta(meta);
    plan.setMachine(machine);
    plan.setNumOfFractions(fractions);
    return plan;
}

ImportedControlPoint controlPoint(double gantry, double weight, std::vector<double> bankA, std::vector<double> bankB) {
    ImportedControlPoint cp;
    cp.gantryAngle = gantry;
    cp.cumulativeWeight = weight;
    cp.bankA = std::move(bankA);
    cp.bankB = std::move(bankB);
    return cp;
}

// Counter-clockwise arc through gantry 0: 10 -> 6 -> 2 -> 358 -> 354, 100 MU per fraction.
ImportedBeam arcThroughZero() {
    ImportedBeam beam;
    beam.name = "Arc";
    beam.isArc = true;
    beam.clockwise = false;
    beam.collimatorAngle = 45.0;
    beam.metersetMU = 100.0;
    beam.isocenter = {1.0, 2.0, 3.0};
    beam.leafBoundaries = {-5.0, 0.0, 5.0};
    const double gantry[] = {10.0, 6.0, 2.0, 358.0, 354.0};
    for (int k = 0; k < 5; ++k)
        beam.controlPoints.push_back(controlPoint(gantry[k], 0.25 * k, {-1.0 * k, -1.0}, {1.0 * k + 1.0, 1.0}));
    return beam;
}

} // namespace

TEST(ImportedPlanConverterTest, ArcSlicesLieBetweenControlPoints) {
    ImportedPlan imported;
    imported.beams.push_back(arcThroughZero());
    const Plan plan = planWith(1000.0, 2);   // 50 machine MU per weight unit

    const VmatSettings settings = ImportedPlanConverter::toVmatSettings(imported);
    ASSERT_EQ(settings.arcs.size(), 1u);
    EXPECT_DOUBLE_EQ(settings.arcs[0].startAngle, 10.0);
    EXPECT_DOUBLE_EQ(settings.arcs[0].stopAngle, 354.0);
    EXPECT_FALSE(settings.arcs[0].clockwise);
    EXPECT_NEAR(settings.arcs[0].span(), 16.0, 1e-9);
    EXPECT_NEAR(settings.controlPointSpacingDeg, 4.0, 1e-9);
    EXPECT_NEAR(settings.optimizationSpacingDeg, 12.0, 1e-9);
    EXPECT_DOUBLE_EQ(settings.arcs[0].collimatorAngle, 45.0);

    const ArcSequenceResult seq = ImportedPlanConverter::toArcSequence(imported, plan);
    ASSERT_EQ(seq.controlPoints.size(), 4u);
    const double expectedOffsets[] = {2.0, 6.0, 10.0, 14.0};
    const double expectedGantry[] = {8.0, 4.0, 0.0, 356.0};
    for (size_t j = 0; j < 4; ++j) {
        const auto& cp = seq.controlPoints[j];
        EXPECT_NEAR(cp.geometry.offset, expectedOffsets[j], 1e-9);
        EXPECT_NEAR(cp.geometry.width, 4.0, 1e-9);
        EXPECT_NEAR(cp.geometry.gantryAngle, expectedGantry[j], 1e-9);
        EXPECT_DOUBLE_EQ(cp.geometry.collimatorAngle, 45.0);
        EXPECT_NEAR(cp.mu, 0.25 * 100.0 / 50.0, 1e-12);
        EXPECT_NEAR(cp.aperture.bankA[0], -(j + 0.5), 1e-12);   // mean of control points j and j+1
        EXPECT_NEAR(cp.aperture.bankB[0], j + 1.5, 1e-12);
        EXPECT_EQ(cp.leafPairBoundariesZ.size(), 3u);
    }
    EXPECT_NEAR(seq.totalMU, 2.0, 1e-9);
    // Two 12-degree target sectors cannot tile 16 degrees, so the whole arc is one sector.
    EXPECT_EQ(seq.controlPoints.back().geometry.sectorIndex, 0u);
}

TEST(ImportedPlanConverterTest, JawsCloseShieldedPairsAndClipTips) {
    ImportedControlPoint cp = controlPoint(0.0, 0.0, {-6.0, -1.0, 0.0}, {6.0, 2.0, 5.0});
    cp.jawX = {-3.0, 4.0};
    cp.jawY = {-5.0, 8.0};
    const Aperture ap = ImportedPlanConverter::collimatedAperture(cp, {-10.0, 0.0, 10.0, 20.0});
    ASSERT_EQ(ap.bankA.size(), 3u);
    EXPECT_DOUBLE_EQ(ap.bankA[0], -3.0);
    EXPECT_DOUBLE_EQ(ap.bankB[0], 4.0);
    EXPECT_DOUBLE_EQ(ap.bankA[1], -1.0);
    EXPECT_DOUBLE_EQ(ap.bankB[1], 2.0);
    EXPECT_DOUBLE_EQ(ap.bankA[2], ap.bankB[2]);   // above the Y2 jaw: closed
}

TEST(ImportedPlanConverterTest, StaticBeamsBecomeStepAndShootSegments) {
    ImportedBeam beam;
    beam.metersetMU = 60.0;
    beam.couchAngle = -10.0;
    beam.collimatorAngle = 45.0;
    beam.leafBoundaries = {-5.0, 5.0};
    beam.controlPoints = {controlPoint(90.0, 0.0, {-2.0}, {2.0}), controlPoint(90.0, 0.5, {-2.0}, {2.0}),
                          controlPoint(90.0, 0.5, {-1.0}, {1.0}), controlPoint(90.0, 1.0, {-1.0}, {1.0})};
    ImportedPlan imported;
    imported.numFractions = 3;
    imported.prescribedDose = 54.0;
    imported.beams = {beam};

    PlanConfig config;
    ImportedPlanConverter::applyToPlanConfig(imported, config);
    EXPECT_EQ(config.deliveryMode, DeliveryMode::SMLC);
    EXPECT_EQ(config.numFractions, 3);
    EXPECT_EQ(config.gantryAngles, std::vector<double>{90.0});
    EXPECT_EQ(config.couchAngles, std::vector<double>{-10.0});
    EXPECT_EQ(config.collimatorAngles, std::vector<double>{45.0});
    ASSERT_TRUE(config.prescribedDose.has_value());
    EXPECT_DOUBLE_EQ(*config.prescribedDose, 54.0);

    const auto seqs = ImportedPlanConverter::toLeafSequences(imported, planWith(1000.0, 3));
    ASSERT_EQ(seqs.size(), 1u);
    ASSERT_EQ(seqs[0].segments.size(), 2u);   // the leaf move between CPs 1 and 2 delivers nothing
    EXPECT_DOUBLE_EQ(seqs[0].segments[1].bankA[0], -1.0);
    // 60 MU/fx x 3 fx over 100 MU per weight unit (SAD 1000)
    EXPECT_NEAR(seqs[0].totalMU, 60.0 * 3.0 / 100.0, 1e-12);
}

TEST(ImportedPlanConverterTest, IsocenterComesFromConfigThenImportedPlanThenTarget) {
    const Plan plan;   // no patient data: computeIsoCenter() yields the origin
    PlanConfig config;
    auto patient = std::make_shared<PatientData>();

    auto iso = PlanBuilder::resolveIsocenter(config, plan, patient.get());
    EXPECT_EQ(iso, (std::array<double, 3>{0.0, 0.0, 0.0}));

    auto imported = std::make_shared<ImportedPlan>();
    imported->beams.push_back(arcThroughZero());
    patient->setImportedPlan(imported);
    iso = PlanBuilder::resolveIsocenter(config, plan, patient.get());
    EXPECT_EQ(iso, (std::array<double, 3>{1.0, 2.0, 3.0}));

    config.isocenter = std::array<double, 3>{4.0, 5.0, 6.0};
    iso = PlanBuilder::resolveIsocenter(config, plan, patient.get());
    EXPECT_EQ(iso, (std::array<double, 3>{4.0, 5.0, 6.0}));
}

TEST(ImportedPlanConverterTest, ArcPlanMismatchesListWhatKeepsAPlanFromDeliveringTheRtplan) {
    ImportedPlan imported;
    imported.numFractions = 2;
    imported.beams.push_back(arcThroughZero());

    Plan plan = planWith(1000.0, 2);
    EXPECT_FALSE(ImportedPlanConverter::arcPlanMismatches(imported, plan).empty());   // not VMAT
    plan.setDeliveryMode(DeliveryMode::VMAT);
    plan.setVmatSettings(ImportedPlanConverter::toVmatSettings(imported));
    plan.stfProperties().isoCenters = {{1.0, 2.0, 3.0}};
    EXPECT_TRUE(ImportedPlanConverter::arcPlanMismatches(imported, plan).empty());

    auto vmat = plan.getVmatSettings();
    vmat.arcs[0].collimatorAngle = 0.0;
    plan.setVmatSettings(vmat);
    plan.setNumOfFractions(3);
    plan.stfProperties().isoCenters = {{1.0, 2.0, 4.0}};
    const auto m = ImportedPlanConverter::arcPlanMismatches(imported, plan);
    ASSERT_EQ(m.size(), 3u);
    EXPECT_NE(m[0].find("collimator"), std::string::npos);
    EXPECT_NE(m[1].find("fraction"), std::string::npos);
    EXPECT_NE(m[2].find("isocenter"), std::string::npos);
}

} // namespace optirad::tests
