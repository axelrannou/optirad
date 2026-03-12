#include <gtest/gtest.h>
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "core/Machine.hpp"
#include "core/Stf.hpp"
#include "core/Beam.hpp"
#include "core/Ray.hpp"
#include "io/MachineLoader.hpp"
#include <cmath>

namespace optirad::tests {

class PhotonIMRTStfGeneratorTest : public ::testing::Test {
protected:
    Machine machine;

    void SetUp() override {
        machine = MachineLoader::load("photons", "Generic");
    }
};

// ============================================================================
// Basic STF generation
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, GenerateStfCreatesCorrectNumberOfBeams) {
    // gantry 0, 40, 80, ..., 320 → 9 beams
    PhotonIMRTStfGenerator gen(0.0, 40.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_EQ(stf.getCount(), 9u);
}

TEST_F(PhotonIMRTStfGeneratorTest, GenerateStfSingleBeam) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_EQ(stf.getCount(), 1u);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectGantryAngles) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 4u); // 0, 90, 180, 270

    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getGantryAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getGantryAngle(), 90.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getGantryAngle(), 180.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getGantryAngle(), 270.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveZeroCouchAngle) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    for (size_t i = 0; i < stf.getCount(); ++i) {
        EXPECT_DOUBLE_EQ(stf.getBeam(i)->getCouchAngle(), 0.0);
    }
}

// ============================================================================
// Beam properties from machine
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectSAD) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getSAD(), 1000.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectSCD) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getSCD(), 500.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectBixelWidth) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 5.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getBixelWidth(), 5.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectIsocenter) {
    Vec3 iso = {-25.7, -289.4, 53.2};
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, iso);
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto& beamIso = stf.getBeam(0)->getIsocenter();
    EXPECT_NEAR(beamIso[0], -25.7, 1e-10);
    EXPECT_NEAR(beamIso[1], -289.4, 1e-10);
    EXPECT_NEAR(beamIso[2], 53.2, 1e-10);
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectRadiationMode) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);
    gen.setRadiationMode("photons");

    Stf stf = gen.generateStf();
    EXPECT_EQ(stf.getBeam(0)->getRadiationMode(), "photons");
}

TEST_F(PhotonIMRTStfGeneratorTest, BeamsHaveCorrectMachineName) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_EQ(stf.getBeam(0)->getMachineName(), "Generic");
}

// ============================================================================
// Rays within beams
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, BeamsContainRays) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    EXPECT_GT(stf.getBeam(0)->getNumOfRays(), 0u);
}

TEST_F(PhotonIMRTStfGeneratorTest, RaysHaveCorrectEnergy) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto* beam = stf.getBeam(0);
    for (size_t i = 0; i < beam->getNumOfRays(); ++i) {
        EXPECT_DOUBLE_EQ(beam->getRay(i)->getEnergy(), 6.0);
    }
}

TEST_F(PhotonIMRTStfGeneratorTest, RaysHaveBeamletCorners) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto* beam = stf.getBeam(0);
    ASSERT_GT(beam->getNumOfRays(), 0u);

    // Check that corners are properly computed (not all zeros for non-center rays)
    const auto* ray = beam->getRay(0);
    const auto& corners = ray->getBeamletCornersAtIso();
    // Corners should form a square of bixelWidth side length
    double dx = corners[0][0] - corners[1][0];
    double dz = corners[1][2] - corners[2][2];
    EXPECT_NEAR(std::abs(dx), 7.0, 1e-10); // bixel width
    EXPECT_NEAR(std::abs(dz), 7.0, 1e-10); // bixel width
}

TEST_F(PhotonIMRTStfGeneratorTest, RaysHaveSCDCorners) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto* beam = stf.getBeam(0);
    const auto* ray = beam->getRay(0);

    // SCD corners should be at y = SCD - SAD = -500 (for gantry=0)
    const auto& scdCorners = ray->getRayCornersSCD();
    for (const auto& c : scdCorners) {
        EXPECT_NEAR(c[1], -500.0, 1e-10);
    }

    // SCD corners should be scaled by SCD/SAD = 0.5 relative to iso corners
    const auto& isoCorners = ray->getBeamletCornersAtIso();
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(scdCorners[i][0], isoCorners[i][0] * 0.5, 1e-10);
        EXPECT_NEAR(scdCorners[i][2], isoCorners[i][2] * 0.5, 1e-10);
    }
}

TEST_F(PhotonIMRTStfGeneratorTest, RaysHaveTargetPoints) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto* beam = stf.getBeam(0);
    const auto* ray = beam->getRay(0);

    // Target point BEV should be (2*x, SAD, 2*z) for a ray at (x, 0, z)
    const auto& pos = ray->getRayPosBev();
    const auto& target = ray->getTargetPointBev();
    EXPECT_NEAR(target[0], 2.0 * pos[0], 1e-10);
    EXPECT_NEAR(target[1], 1000.0, 1e-10);
    EXPECT_NEAR(target[2], 2.0 * pos[2], 1e-10);
}

// ============================================================================
// Source point
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, SourcePointBev) {
    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    const auto& sp = stf.getBeam(0)->getSourcePointBev();
    EXPECT_DOUBLE_EQ(sp[0], 0.0);
    EXPECT_DOUBLE_EQ(sp[1], -1000.0);
    EXPECT_DOUBLE_EQ(sp[2], 0.0);
}

// ============================================================================
// Totals across all beams
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, TotalRaysAcrossBeams) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    size_t totalRays = stf.getTotalNumOfRays();
    
    // All beams should have the same number of rays (same field size)
    size_t raysPerBeam = stf.getBeam(0)->getNumOfRays();
    EXPECT_EQ(totalRays, raysPerBeam * 4);  // 4 beams
}

TEST_F(PhotonIMRTStfGeneratorTest, TotalBixelsAcrossBeams) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    Stf stf = gen.generateStf();
    // For photons: total bixels = total rays (1 bixel per ray)
    EXPECT_EQ(stf.getTotalNumOfBixels(), stf.getTotalNumOfRays());
}

// ============================================================================
// Legacy StfProperties generation still works
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, LegacyGenerateStillWorks) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {1.0, 2.0, 3.0});
    auto props = gen.generate();

    ASSERT_NE(props, nullptr);
    EXPECT_EQ(props->numOfBeams, 4u);
    EXPECT_DOUBLE_EQ(props->bixelWidth, 7.0);
    EXPECT_EQ(props->gantryAngles.size(), 4u);
    EXPECT_EQ(props->isoCenters.size(), 4u);
}

// ============================================================================
// Target-aware STF generation
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, TargetAwareGenerationProducesDifferentRayCounts) {
    // Create a target elongated along the x-axis. 
    // At different gantry angles, the projection of this target onto the 
    // BEV isocenter plane will have different extents, producing different ray counts.
    Vec3 iso = {0.0, 0.0, 0.0};

    // Target voxels spread along x: from -35 to +35 mm, narrow in y and z
    std::vector<Vec3> targetCoords;
    for (double x = -35.0; x <= 35.0; x += 7.0) {
        targetCoords.push_back({x, 0.0, 0.0});
    }

    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, iso);
    gen.setMachine(machine);
    gen.setTargetVoxelCoords(targetCoords);

    Stf stf = gen.generateStf();

    ASSERT_EQ(stf.getCount(), 4u); // 0, 90, 180, 270

    // Beams at 0° and 180° see the target spread in BEV-x → more rays
    // Beams at 90° and 270° see the target spread along BEV-y (depth) → fewer rays
    size_t rays0 = stf.getBeam(0)->getNumOfRays();
    size_t rays90 = stf.getBeam(1)->getNumOfRays();
    size_t rays180 = stf.getBeam(2)->getNumOfRays();
    size_t rays270 = stf.getBeam(3)->getNumOfRays();

    // 0° and 180° should have more rays than 90° and 270°
    EXPECT_GT(rays0, rays90);
    EXPECT_GT(rays180, rays270);

    // The ray counts should NOT all be identical
    bool allSame = (rays0 == rays90) && (rays90 == rays180) && (rays180 == rays270);
    EXPECT_FALSE(allSame);
}

TEST_F(PhotonIMRTStfGeneratorTest, TargetAwareBeamsStillHaveCorners) {
    Vec3 iso = {0.0, 0.0, 0.0};
    std::vector<Vec3> targetCoords;
    for (double x = -14.0; x <= 14.0; x += 7.0) {
        for (double z = -14.0; z <= 14.0; z += 7.0) {
            targetCoords.push_back({x, 0.0, z});
        }
    }

    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 7.0, iso);
    gen.setMachine(machine);
    gen.setTargetVoxelCoords(targetCoords);

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 1u);

    const auto* beam = stf.getBeam(0);
    ASSERT_GT(beam->getNumOfRays(), 0u);

    // Every ray should have well-formed corners
    for (size_t i = 0; i < beam->getNumOfRays(); ++i) {
        const auto* ray = beam->getRay(i);
        ASSERT_NE(ray, nullptr);
        const auto& corners = ray->getBeamletCornersAtIso();
        // The corners of each beamlet should span bixelWidth in x and z
        double dx = corners[0][0] - corners[1][0]; // (+x) - (-x) = bixelWidth
        double dz = corners[1][2] - corners[2][2]; // (+z) - (-z) = bixelWidth
        EXPECT_NEAR(std::abs(dx), 7.0, 1e-10);
        EXPECT_NEAR(std::abs(dz), 7.0, 1e-10);
    }
}

TEST_F(PhotonIMRTStfGeneratorTest, TargetAwareWithCTResolution) {
    Vec3 iso = {0.0, 0.0, 0.0};
    std::vector<Vec3> targetCoords = {{0.0, 0.0, 0.0}};

    PhotonIMRTStfGenerator gen(0.0, 360.0, 360.0, 3.0, iso);
    gen.setMachine(machine);
    gen.setTargetVoxelCoords(targetCoords);
    gen.setCTResolution({7.0, 7.0, 7.0});

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 1u);

    // With bixelWidth=3, no Grid/StructureSet → no 3D dilation, single voxel produces 1 ray.
    // Then ctRes padding: 3 < 7, pad = floor(7/3) = 2 → 5x5 = 25 rays
    EXPECT_EQ(stf.getBeam(0)->getNumOfRays(), 25u);
}

TEST_F(PhotonIMRTStfGeneratorTest, FallbackToGridWhenNoTargetCoords) {
    // When no target voxel coords are set, should fall back to grid-based generation
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);
    // Don't call setTargetVoxelCoords

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 4u);

    // All beams should have the same ray count (grid-based, same field size)
    size_t rays0 = stf.getBeam(0)->getNumOfRays();
    for (size_t i = 1; i < stf.getCount(); ++i) {
        EXPECT_EQ(stf.getBeam(i)->getNumOfRays(), rays0);
    }
}

// ============================================================================
// Couch angle support
// ============================================================================

TEST_F(PhotonIMRTStfGeneratorTest, CouchStartStepStop) {
    // 4 gantry angles: 0, 90, 180, 270
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    // 2 couch angles: 0, 15 → Cartesian product: 4 gantry × 2 couch = 8 beams
    gen.setCouchAngles(0.0, 15.0, 30.0);

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 8u);

    // Arc 1 (couch=0): gantry 0, 90, 180, 270
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getGantryAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getCouchAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getGantryAngle(), 90.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getCouchAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getGantryAngle(), 180.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getCouchAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getGantryAngle(), 270.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getCouchAngle(), 0.0);

    // Arc 2 (couch=15): gantry 0, 90, 180, 270
    EXPECT_DOUBLE_EQ(stf.getBeam(4)->getGantryAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(4)->getCouchAngle(), 15.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(5)->getGantryAngle(), 90.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(5)->getCouchAngle(), 15.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(6)->getGantryAngle(), 180.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(6)->getCouchAngle(), 15.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(7)->getGantryAngle(), 270.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(7)->getCouchAngle(), 15.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, CouchStartStepStopProps) {
    // Also verify the generate() path (StfProperties) does Cartesian product
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);
    gen.setCouchAngles(0.0, 10.0, 30.0);  // 3 couch: 0, 10, 20 → 4×3 = 12 beams

    auto props = gen.generate();
    ASSERT_EQ(props->numOfBeams, 12u);
    ASSERT_EQ(props->gantryAngles.size(), 12u);
    ASSERT_EQ(props->couchAngles.size(), 12u);
    EXPECT_TRUE(props->isValid());

    // First arc: gantry=[0,90,180,270] couch=0
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(props->couchAngles[i], 0.0);
    }
    // Second arc: gantry=[0,90,180,270] couch=10
    for (size_t i = 4; i < 8; ++i) {
        EXPECT_DOUBLE_EQ(props->couchAngles[i], 10.0);
    }
    // Third arc: gantry=[0,90,180,270] couch=20
    for (size_t i = 8; i < 12; ++i) {
        EXPECT_DOUBLE_EQ(props->couchAngles[i], 20.0);
    }
}

TEST_F(PhotonIMRTStfGeneratorTest, CouchAnglesPairedWithGantry) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    // Explicit couch angle list (must be same length as gantry)
    gen.setCouchAngles({10.0, 20.0, 30.0, 40.0});

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 4u);

    // Verify paired: beam[i] = (gantry[i], couch[i])
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getGantryAngle(), 0.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getCouchAngle(), 10.0);

    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getGantryAngle(), 90.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getCouchAngle(), 20.0);

    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getGantryAngle(), 180.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getCouchAngle(), 30.0);

    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getGantryAngle(), 270.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getCouchAngle(), 40.0);
}

TEST_F(PhotonIMRTStfGeneratorTest, UniformCouchAngle) {
    // When couch step=0, all beams get the same couch angle
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    // Single couch angle = 25.0 for all beams
    gen.setCouchAngles(25.0, 0.0, 0.0);

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 4u);

    for (size_t i = 0; i < stf.getCount(); ++i) {
        EXPECT_DOUBLE_EQ(stf.getBeam(i)->getCouchAngle(), 25.0);
    }
}

TEST_F(PhotonIMRTStfGeneratorTest, GeneratePropsWithCouchAngles) {
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);
    gen.setCouchAngles({5.0, 10.0, 15.0, 20.0});

    auto props = gen.generate();
    ASSERT_EQ(props->gantryAngles.size(), 4u);
    ASSERT_EQ(props->couchAngles.size(), 4u);

    EXPECT_DOUBLE_EQ(props->couchAngles[0], 5.0);
    EXPECT_DOUBLE_EQ(props->couchAngles[1], 10.0);
    EXPECT_DOUBLE_EQ(props->couchAngles[2], 15.0);
    EXPECT_DOUBLE_EQ(props->couchAngles[3], 20.0);
    EXPECT_TRUE(props->isValid());
}

TEST_F(PhotonIMRTStfGeneratorTest, CouchAnglesResizedToMatchGantry) {
    // If couch list is shorter, it gets padded with last value
    PhotonIMRTStfGenerator gen(0.0, 90.0, 360.0, 7.0, {0.0, 0.0, 0.0});
    gen.setMachine(machine);

    // Only 2 couch angles for 4 gantry angles
    gen.setCouchAngles({10.0, 20.0});

    Stf stf = gen.generateStf();
    ASSERT_EQ(stf.getCount(), 4u);

    EXPECT_DOUBLE_EQ(stf.getBeam(0)->getCouchAngle(), 10.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(1)->getCouchAngle(), 20.0);
    // Padded with last value (20.0)
    EXPECT_DOUBLE_EQ(stf.getBeam(2)->getCouchAngle(), 20.0);
    EXPECT_DOUBLE_EQ(stf.getBeam(3)->getCouchAngle(), 20.0);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
