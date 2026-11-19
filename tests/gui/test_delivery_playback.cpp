#include <gtest/gtest.h>
#include "gui/AppState.hpp"
#include <cmath>
#include <limits>

namespace optirad::tests {
namespace {

class DeliveryPlaybackTest : public testing::Test {
protected:
    GuiAppState state;
    std::shared_ptr<ArcPipelineResult> result;

    void SetUp() override {
        state.plan = std::make_shared<Plan>();
        state.plan->setDeliveryMode(DeliveryMode::VMAT);
        VmatSettings settings;
        settings.arcs = {{350.0, 10.0, true, 0.0}, {10.0, 350.0, false, 30.0}};
        state.plan->setVmatSettings(settings);
        state.stf = std::make_shared<Stf>();
        result = std::make_shared<ArcPipelineResult>();
        result->controlPointStf = std::make_shared<Stf>();
        const double angles[] = {355.0, 5.0, 5.0, 355.0};
        const double seconds[] = {1.0, 2.0, 1.0, 3.0};
        for (size_t i = 0; i < 4; ++i) {
            Beam beam;
            beam.setIsocenter(15.0, 20.0, 25.0);
            beam.setGantryAngle(angles[i]);
            beam.setCouchAngle(i < 2 ? 0.0 : 30.0);
            beam.setSAD(1000.0);
            beam.setSCD(500.0);
            result->controlPointStf->addBeam(beam);
            DeliveredControlPoint cp;
            cp.geometry.arcIndex = i / 2;
            cp.geometry.sectorIndex = i;
            cp.geometry.gantryAngle = angles[i];
            cp.geometry.couchAngle = beam.getCouchAngle();
            cp.geometry.offset = (i % 2 + 0.5) * 10.0;
            cp.geometry.width = 10.0;
            cp.timeSec = seconds[i];
            cp.mu = 1.0;
            cp.gantrySpeed = 10.0 / seconds[i];
            cp.doseRate = 1.0;
            cp.aperture.bankA = {-20.0 - i * 10.0, -10.0};
            cp.aperture.bankB = {20.0 + i * 10.0, 10.0};
            cp.leafPairBoundariesZ = {-5.0, 0.0, 5.0};
            result->sequence.controlPoints.push_back(cp);
            result->sequence.totalTimeSec += cp.timeSec;
            state.stf->addBeam(beam);
        }
        const int id = state.doseStore.addDose("Delivered arcs", nullptr, nullptr);
        state.arcCache[id].result = result;
        state.delivery.synchronize(state);
    }
};

TEST_F(DeliveryPlaybackTest, WholePlanClockCrossesArcsAndStopsAtTheEnd) {
    auto& delivery = state.delivery;
    ASSERT_TRUE(delivery.available());
    EXPECT_EQ(delivery.groupCount(), 2u);
    EXPECT_DOUBLE_EQ(delivery.duration(), 7.0);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 350.0);
    delivery.togglePlaying();
    delivery.advance(0.5);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 355.0);
    delivery.advance(2.5);
    EXPECT_EQ(delivery.frame().group, 1u);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 10.0);
    delivery.advance(0.5);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 5.0);
    delivery.playbackSpeed = 2.0f;
    delivery.advance(20.0);
    EXPECT_DOUBLE_EQ(delivery.time(), delivery.duration());
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 350.0);
    EXPECT_FALSE(delivery.playing);
    delivery.togglePlaying();
    EXPECT_DOUBLE_EQ(delivery.time(), 0.0);
    EXPECT_TRUE(delivery.playing);
    delivery.rewind();
    EXPECT_FALSE(delivery.playing);
}

TEST_F(DeliveryPlaybackTest, ScrubbingIsContinuousAtWrapAndLeafBoundariesButNotAcrossArcs) {
    auto& delivery = state.delivery;
    delivery.seek(1.0 - 1e-8);
    const auto before = delivery.frame().aperture;
    EXPECT_NEAR(std::remainder(delivery.frame().beam.getGantryAngle(), 360.0), 0.0, 1e-6);
    delivery.seek(1.0 + 1e-8);
    EXPECT_NEAR(delivery.frame().aperture.bankA[0], before.bankA[0], 1e-6);
    delivery.seek(2.0);
    EXPECT_DOUBLE_EQ(delivery.frame().aperture.bankA[0], -30.0);
    delivery.selectGroup(1);
    EXPECT_DOUBLE_EQ(delivery.time(), 3.0);
    EXPECT_DOUBLE_EQ(delivery.frame().aperture.bankA[0], -40.0);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getCouchAngle(), 30.0);
    EXPECT_DOUBLE_EQ(delivery.frame().doseRate, state.plan->machineMuPerWeight());
}

TEST_F(DeliveryPlaybackTest, ArcAndLeafGeometryUseTheSameLpsSourceAndCollimatorPlane) {
    auto& delivery = state.delivery;
    delivery.seek(3.5);
    const auto& beam = delivery.frame().beam;
    const auto inverse = transpose(getRotationMatrix(beam.getGantryAngle(), beam.getCouchAngle()));
    ASSERT_FALSE(delivery.overlay().leaves.empty());
    const auto tip = inverse * vecSub(delivery.overlay().leaves.front().to, beam.getIsocenter());
    EXPECT_NEAR(tip[0], delivery.frame().aperture.bankA[0] * beam.getSCD() / beam.getSAD(), 1e-9);
    for (const auto& line : delivery.overlay().leaves) {
        for (const auto& p : {line.from, line.to}) {
            const auto local = inverse * vecSub(p, beam.getIsocenter());
            EXPECT_NEAR(local[1], beam.getSCD() - beam.getSAD(), 1e-9);
        }
    }
    const auto source = vecAdd(beam.getIsocenter(),
        getRotationMatrix(beam.getGantryAngle(), beam.getCouchAngle()) * Vec3{0.0, -beam.getSAD(), 0.0});
    const auto& ray = delivery.overlay().aim.front();
    EXPECT_LT(norm(cross(vecSub(ray.from, source), vecSub(ray.to, source))), 1e-7);
    ASSERT_GT(delivery.overlay().path.size(), 4u);
    EXPECT_NEAR(norm(vecSub(delivery.overlay().path.front().from, beam.getIsocenter())), 500.0, 1e-8);
    // Each 20-degree arc: 10 two-degree segments, 2 arrow lines, then one tick per control point.
    for (size_t arrow : {10u, 24u}) {
        const auto& end = delivery.overlay().path[arrow - 1];
        const auto& head = delivery.overlay().path[arrow];
        EXPECT_EQ(end.to, head.from);
        EXPECT_LT(dot(vecSub(head.to, head.from), vecSub(end.to, end.from)), 0.0);
    }
    const auto fit = delivery.overlay().fitPoints;
    delivery.seek(6.5);
    EXPECT_EQ(delivery.overlay().fitPoints, fit);
}

TEST_F(DeliveryPlaybackTest, VmatArcShowsOneRadialTickPerControlPoint) {
    const auto& path = state.delivery.overlay().path;
    ASSERT_EQ(path.size(), 28u);
    const auto& points = result->sequence.controlPoints;
    const size_t tickIndex[] = {12u, 13u, 26u, 27u};
    for (size_t i = 0; i < 4; ++i) {
        const auto& tick = path[tickIndex[i]];
        const auto beam = *result->controlPointStf->getBeam(i);
        const auto iso = beam.getIsocenter();
        const auto mid = vecScale(vecAdd(tick.from, tick.to), 0.5);
        EXPECT_NEAR(norm(vecSub(mid, iso)), 500.0, 1e-8);
        EXPECT_NEAR(norm(vecSub(tick.to, tick.from)), 16.0, 1e-8);
        EXPECT_LT(norm(cross(vecSub(tick.to, tick.from), vecSub(mid, iso))), 1e-6);
        const auto expected = vecAdd(iso, getRotationMatrix(points[i].geometry.gantryAngle,
            points[i].geometry.couchAngle) * Vec3{0.0, -500.0, 0.0});
        EXPECT_LT(norm(vecSub(mid, expected)), 1e-8);
    }
}

TEST_F(DeliveryPlaybackTest, InPlaceRetimingIsDetectedEvenWhenTotalDurationDoesNotChange) {
    result->sequence.controlPoints[0].timeSec = 2.0;
    result->sequence.controlPoints[1].timeSec = 1.0;
    state.delivery.synchronize(state);
    EXPECT_DOUBLE_EQ(state.delivery.duration(), 7.0);
    EXPECT_DOUBLE_EQ(state.delivery.frame().gantrySpeed, 5.0);
    state.delivery.seek(1.0);
    EXPECT_DOUBLE_EQ(state.delivery.frame().beam.getGantryAngle(), 355.0);
}

TEST_F(DeliveryPlaybackTest, GeometryPreviewWorksBeforeMlcAndRespectsFullArcDirection) {
    state.doseStore.clear();
    VmatSettings settings;
    settings.arcs = {{0.0, 360.0, false, 25.0}};
    state.plan->setVmatSettings(settings);
    state.delivery.synchronize(state);
    auto& delivery = state.delivery;
    ASSERT_TRUE(delivery.available());
    EXPECT_TRUE(delivery.estimated());
    EXPECT_FALSE(delivery.frame().hasAperture);
    EXPECT_DOUBLE_EQ(delivery.duration(), 60.0);
    delivery.seek(15.0);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 270.0);
    EXPECT_FALSE(delivery.overlay().path.empty());
    EXPECT_EQ(delivery.overlay().leaves.size(), 4u);
}

TEST_F(DeliveryPlaybackTest, SmlcUsesSelectedCachedSegmentsWithBeamOffTransitionsAndCalibratedTiming) {
    state.doseStore.clear();
    state.plan->setDeliveryMode(DeliveryMode::SMLC);
    state.plan->setNumOfFractions(1);
    state.stf = std::make_shared<Stf>();
    Beam beam;
    beam.setGantryAngle(0.0);
    state.stf->addBeam(beam);
    beam.setGantryAngle(90.0);
    state.stf->addBeam(beam);
    LeafSequenceResult sequence;
    Aperture first, second;
    first.bankA = {-20.0}; first.bankB = {20.0}; first.weight = 0.1;
    second.bankA = {-10.0}; second.bankB = {10.0}; second.weight = 0.2;
    sequence.segments = {first, second};
    sequence.leafPairBoundariesZ = {-5.0, 5.0};
    const int id = state.doseStore.addDose("Sequenced", nullptr, nullptr);
    state.seqCache[id].sequences = {sequence, sequence};
    state.leafSequences.clear(); // playback must use the selected cache, not unrelated mutable GUI data
    state.delivery.synchronize(state);
    auto& delivery = state.delivery;
    ASSERT_TRUE(delivery.available());
    EXPECT_FALSE(delivery.isVmat());
    EXPECT_TRUE(delivery.estimated());
    EXPECT_EQ(delivery.frame().segment, 0);
    delivery.seek(1.0 + 1.0 / 12.0);
    EXPECT_FALSE(delivery.frame().beamOn);
    EXPECT_NEAR(delivery.frame().aperture.bankA[0], -15.0, 1e-8);
    delivery.selectSegment(0, 1);
    EXPECT_TRUE(delivery.frame().beamOn);
    EXPECT_EQ(delivery.frame().segment, 1);
    delivery.selectGroup(1);
    EXPECT_EQ(delivery.frame().beamIndex, 1u);
    EXPECT_FALSE(delivery.frame().beamOn);
    delivery.selectSegment(1, 0);
    EXPECT_TRUE(delivery.frame().beamOn);
    EXPECT_DOUBLE_EQ(delivery.frame().beam.getGantryAngle(), 90.0);
    EXPECT_DOUBLE_EQ(delivery.frame().doseRate, 10.0);
}

TEST_F(DeliveryPlaybackTest, SourceReplacementResetsPlaybackWithoutResettingVisibilityOrSpeed) {
    auto& delivery = state.delivery;
    delivery.showOverlays = false;
    delivery.playbackSpeed = 2.5f;
    delivery.togglePlaying();
    delivery.advance(1.0);
    const double time = delivery.time();
    delivery.synchronize(state);
    EXPECT_DOUBLE_EQ(delivery.time(), time);
    EXPECT_TRUE(delivery.playing);
    state.patientData = std::make_shared<PatientData>();
    delivery.synchronize(state);
    EXPECT_DOUBLE_EQ(delivery.time(), 0.0);
    EXPECT_FALSE(delivery.playing);
    EXPECT_FALSE(delivery.showOverlays);
    EXPECT_FLOAT_EQ(delivery.playbackSpeed, 2.5f);
    state.plan.reset();
    state.stf.reset();
    state.doseStore.clear();
    delivery.synchronize(state);
    EXPECT_FALSE(delivery.available());
    EXPECT_FALSE(delivery.frame().valid);
    EXPECT_TRUE(delivery.overlay().leaves.empty());
}

TEST_F(DeliveryPlaybackTest, MalformedMlcIsReportedAndClearsStaleOverlay) {
    result->sequence.controlPoints[0].leafPairBoundariesZ.clear();
    state.doseStore.renameDose(0, "Invalid MLC");
    state.delivery.synchronize(state);
    EXPECT_FALSE(state.delivery.available());
    EXPECT_FALSE(state.delivery.frame().valid);
    EXPECT_TRUE(state.delivery.overlay().leaves.empty());
    EXPECT_NE(state.delivery.description().find("inconsistent MLC"), std::string::npos);
    EXPECT_THROW(state.delivery.seek(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    EXPECT_THROW(state.delivery.advance(-1.0), std::invalid_argument);
}

TEST_F(DeliveryPlaybackTest, NonfiniteLeafCoordinatesAreReportedWithoutBuildingInvalidLines) {
    result->sequence.controlPoints[1].aperture.bankB[0] = std::numeric_limits<double>::quiet_NaN();
    state.doseStore.renameDose(0, "Invalid coordinates");
    state.delivery.synchronize(state);
    EXPECT_FALSE(state.delivery.available());
    EXPECT_TRUE(state.delivery.overlay().leaves.empty());
    EXPECT_NE(state.delivery.description().find("invalid MLC coordinates"), std::string::npos);
}

TEST_F(DeliveryPlaybackTest, ZeroDurationDoesNotBypassGeometryValidation) {
    result->sequence.controlPoints[0].timeSec = 0.0;
    result->sequence.controlPoints[0].geometry.width = std::numeric_limits<double>::quiet_NaN();
    state.delivery.synchronize(state);
    EXPECT_FALSE(state.delivery.available());
    EXPECT_TRUE(state.delivery.overlay().path.empty());
    EXPECT_NE(state.delivery.description().find("invalid beam coordinates"), std::string::npos);
}

} // namespace
} // namespace optirad::tests
