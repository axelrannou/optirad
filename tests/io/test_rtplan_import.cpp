#include <gtest/gtest.h>
#include "io/DicomExporter.hpp"
#include "io/DicomImporter.hpp"
#include "io/MachineLoader.hpp"
#include "core/workflow/ImportedPlanConverter.hpp"
#include "core/Beam.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"

#include <filesystem>

namespace optirad::tests {

TEST(RTPlanImportTest, ExportedStepAndShootPlanRoundTrips) {
    Plan plan;
    Machine machine;
    MachineMeta meta;
    meta.name = "Test";
    meta.SAD = 1000.0;
    machine.setMeta(meta);
    plan.setMachine(machine);
    plan.setNumOfFractions(30);

    Stf stf;
    Beam beam;
    beam.setGantryAngle(40.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({-12.5, 30.0, 7.25});
    stf.addBeam(beam);

    LeafSequenceResult seq;
    seq.leafPairBoundariesZ = {-5.0, 0.0, 5.0};
    const double weights[] = {10.0, 20.0};
    for (int s = 0; s < 2; ++s) {
        Aperture segment;
        segment.bankA = {-10.0 + s, -8.0};
        segment.bankB = {10.0, 6.0 - s};
        segment.weight = weights[s];
        segment.segmentIndex = s;
        seq.segments.push_back(segment);
    }
    seq.totalMU = 30.0;

    const auto dir = std::filesystem::temp_directory_path() / "optirad_test_rtplan_roundtrip";
    std::filesystem::remove_all(dir);
    std::string uid;
    ASSERT_TRUE(DicomExporter().exportRTPlan(plan, stf, {seq}, DicomContext{}, dir.string(), &uid));

    DicomImporter importer;
    const auto imported = importer.importRTPlan((dir / ("RP." + uid + ".dcm")).string());
    ASSERT_NE(imported, nullptr);
    EXPECT_EQ(imported->numFractions, 30);
    ASSERT_EQ(imported->beams.size(), 1u);
    const ImportedBeam& b = imported->beams[0];
    EXPECT_FALSE(b.isArc);
    EXPECT_NEAR(b.controlPoints.front().gantryAngle, 40.0, 1e-6);
    EXPECT_NEAR(b.isocenter[0], -12.5, 1e-3);
    EXPECT_NEAR(b.isocenter[2], 7.25, 1e-3);
    EXPECT_NEAR(b.controlPoints.back().cumulativeWeight, 1.0, 1e-9);
    EXPECT_EQ(b.leafBoundaries, seq.leafPairBoundariesZ);

    const auto back = ImportedPlanConverter::toLeafSequences(*imported, plan);
    ASSERT_EQ(back.size(), 1u);
    ASSERT_EQ(back[0].segments.size(), 2u);
    for (int s = 0; s < 2; ++s) {
        EXPECT_NEAR(back[0].segments[s].weight, weights[s], 1e-3);
        EXPECT_NEAR(back[0].segments[s].bankA[0], seq.segments[s].bankA[0], 1e-6);
        EXPECT_NEAR(back[0].segments[s].bankB[1], seq.segments[s].bankB[1], 1e-6);
    }
    std::filesystem::remove_all(dir);
}

TEST(RTPlanImportTest, ClinicalVmatLungPlan) {
    const std::filesystem::path dir = std::filesystem::path(OPTIRAD_DATA_DIR) / "dicom" / "JOHN_DOE_LUNG";
    if (!std::filesystem::exists(dir)) GTEST_SKIP() << "JOHN_DOE_LUNG test data not available";

    DicomImporter importer;
    ASSERT_TRUE(importer.loadDirectory(dir.string()));
    const auto imported = importer.importRTPlan();
    ASSERT_NE(imported, nullptr);

    EXPECT_EQ(imported->numFractions, 33);
    EXPECT_EQ(imported->patientPosition, "HFS");
    ASSERT_EQ(imported->beams.size(), 2u);
    EXPECT_FALSE(imported->hasMultipleIsocenters());
    for (const auto& b : imported->beams) {
        EXPECT_TRUE(b.isArc);
        EXPECT_FALSE(b.fff);
        EXPECT_EQ(b.leafBoundaries.size(), 61u);
        EXPECT_EQ(b.controlPoints.front().bankA.size(), 60u);
        EXPECT_NEAR(b.controlPoints.back().cumulativeWeight, 1.0, 1e-9);
        EXPECT_NEAR(b.collimatorAngle, 45.0, 1e-6);
    }
    EXPECT_NEAR(imported->beams[0].metersetMU, 206.0, 1.0);
    EXPECT_NEAR(imported->beams[1].metersetMU, 204.5, 1.0);

    const VmatSettings vmat = ImportedPlanConverter::toVmatSettings(*imported);
    ASSERT_EQ(vmat.arcs.size(), 2u);
    EXPECT_NEAR(vmat.arcs[0].span(), 208.0, 0.5);
    EXPECT_NE(vmat.arcs[0].clockwise, vmat.arcs[1].clockwise);

    Plan plan;
    plan.setMachine(MachineLoader::load("photons", "HD120"));
    plan.setNumOfFractions(imported->numFractions);
    const ArcSequenceResult seq = ImportedPlanConverter::toArcSequence(*imported, plan);
    ASSERT_FALSE(seq.controlPoints.empty());
    EXPECT_DOUBLE_EQ(seq.controlPoints.front().geometry.collimatorAngle, 45.0);
    const double expectedWeight = (206.0 + 204.5) / plan.machineMuPerWeight();
    EXPECT_NEAR(seq.totalMU, expectedWeight, 0.01 * expectedWeight);
    EXPECT_GT(seq.totalTimeSec, 0.0);

    // The machine file matches the plan's MLC; its 45-degree collimator angle is supported.
    for (const auto& w : ImportedPlanConverter::validate(*imported, plan.getMachine())) {
        EXPECT_EQ(w.find("leaf"), std::string::npos) << w;
        EXPECT_EQ(w.find("SAD"), std::string::npos) << w;
        EXPECT_EQ(w.find("collimator"), std::string::npos) << w;
    }
}

} // namespace optirad::tests
