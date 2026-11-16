#include <gtest/gtest.h>
#include "io/DicomExporter.hpp"
#include "core/Beam.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"

#include <dcmtk/dcmdata/dctk.h>
#include <filesystem>

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

} // namespace

TEST(DicomExporterTest, MachineMuPerWeightFollowsTheMatRadCalibration) {
    EXPECT_NEAR(planWith(1000.0, 30).machineMuPerWeight(), 100.0 / 30.0, 1e-12);
    EXPECT_NEAR(planWith(800.0, 1).machineMuPerWeight(), 64.0, 1e-12);
    EXPECT_NEAR(planWith(0.0, 0).machineMuPerWeight(), 100.0, 1e-12);   // unknown SAD / fractions
}

TEST(DicomExporterTest, BeamMetersetIsMachineMuPerFraction) {
    const Plan plan = planWith(1000.0, 30);
    Stf stf;
    Beam beam;
    beam.setGantryAngle(0.0);
    beam.setSAD(1000.0);
    stf.addBeam(beam);

    LeafSequenceResult seq;
    seq.beamIndex = 0;
    for (int s = 0; s < 2; ++s) {
        Aperture segment;
        segment.bankA = {-10.0, -10.0};
        segment.bankB = {10.0, 10.0};
        segment.weight = 15.0;
        segment.segmentIndex = s;
        seq.segments.push_back(segment);
    }
    seq.totalMU = 30.0;   // weight units: whole course
    seq.leafPairBoundariesZ = {-5.0, 0.0, 5.0};

    const auto dir = std::filesystem::temp_directory_path() / "optirad_test_rtplan_export";
    std::filesystem::remove_all(dir);
    std::string uid;
    ASSERT_TRUE(DicomExporter().exportRTPlan(plan, stf, {seq}, DicomContext{}, dir.string(), &uid));
    const auto file = dir / ("RP." + uid + ".dcm");
    ASSERT_TRUE(std::filesystem::exists(file));

    DcmFileFormat format;
    ASSERT_TRUE(format.loadFile(file.c_str()).good());
    DcmItem* fractionGroup = nullptr;
    DcmItem* referencedBeam = nullptr;
    auto* dataset = format.getDataset();
    ASSERT_TRUE(dataset->findAndGetSequenceItem(DCM_FractionGroupSequence, fractionGroup, 0).good());
    ASSERT_TRUE(fractionGroup->findAndGetSequenceItem(DCM_ReferencedBeamSequence, referencedBeam, 0).good());
    OFString meterset;
    ASSERT_TRUE(referencedBeam->findAndGetOFString(DCM_BeamMeterset, meterset).good());
    EXPECT_NEAR(std::stod(meterset.c_str()), 100.0, 1e-3);   // 30 x 100 MU / 30 fractions

    // Control point weights stay relative: the last cumulative weight equals the final one.
    DcmItem* beamItem = nullptr;
    ASSERT_TRUE(dataset->findAndGetSequenceItem(DCM_BeamSequence, beamItem, 0).good());
    OFString finalWeight, lastWeight;
    ASSERT_TRUE(beamItem->findAndGetOFString(DCM_FinalCumulativeMetersetWeight, finalWeight).good());
    DcmItem* lastControlPoint = nullptr;
    ASSERT_TRUE(beamItem->findAndGetSequenceItem(DCM_ControlPointSequence, lastControlPoint, -1).good());
    ASSERT_TRUE(lastControlPoint->findAndGetOFString(DCM_CumulativeMetersetWeight, lastWeight).good());
    EXPECT_NEAR(std::stod(lastWeight.c_str()), std::stod(finalWeight.c_str()), 1e-6);
    std::filesystem::remove_all(dir);
}

} // namespace optirad::tests
