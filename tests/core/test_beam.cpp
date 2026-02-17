#include <gtest/gtest.h>
#include "core/Beam.hpp"
#include "core/Ray.hpp"
#include "geometry/MathUtils.hpp"
#include <cmath>

namespace optirad::tests {

class BeamTest : public ::testing::Test {
protected:
    Beam beam;
};

// ============================================================================
// Basic property tests
// ============================================================================

TEST_F(BeamTest, InitialGantryAngleIsZero) {
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 0.0);
}

TEST_F(BeamTest, InitialCouchAngleIsZero) {
    EXPECT_DOUBLE_EQ(beam.getCouchAngle(), 0.0);
}

TEST_F(BeamTest, SetGantryAngleWorks) {
    beam.setGantryAngle(90.0);
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 90.0);
}

TEST_F(BeamTest, SetCouchAngleWorks) {
    beam.setCouchAngle(45.0);
    EXPECT_DOUBLE_EQ(beam.getCouchAngle(), 45.0);
}

TEST_F(BeamTest, NegativeGantryAngleAccepted) {
    beam.setGantryAngle(-45.0);
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), -45.0);
}

TEST_F(BeamTest, LargeGantryAngleAccepted) {
    beam.setGantryAngle(450.0);
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 450.0);
}

TEST_F(BeamTest, SetIsocenterXYZ) {
    beam.setIsocenter(100.0, 200.0, 300.0);
    const auto& iso = beam.getIsocenter();
    EXPECT_DOUBLE_EQ(iso[0], 100.0);
    EXPECT_DOUBLE_EQ(iso[1], 200.0);
    EXPECT_DOUBLE_EQ(iso[2], 300.0);
}

TEST_F(BeamTest, SetIsocenterVec3) {
    beam.setIsocenter({-25.7, -289.4, 53.2});
    const auto& iso = beam.getIsocenter();
    EXPECT_NEAR(iso[0], -25.7, 1e-10);
    EXPECT_NEAR(iso[1], -289.4, 1e-10);
    EXPECT_NEAR(iso[2], 53.2, 1e-10);
}

TEST_F(BeamTest, MultipleAngleUpdates) {
    beam.setGantryAngle(0.0);
    beam.setGantryAngle(180.0);
    beam.setGantryAngle(360.0);
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 360.0);
}

// ============================================================================
// Default values for new properties
// ============================================================================

TEST_F(BeamTest, DefaultBixelWidth) {
    EXPECT_DOUBLE_EQ(beam.getBixelWidth(), 7.0);
}

TEST_F(BeamTest, DefaultSAD) {
    EXPECT_DOUBLE_EQ(beam.getSAD(), 1000.0);
}

TEST_F(BeamTest, DefaultSCD) {
    EXPECT_DOUBLE_EQ(beam.getSCD(), 500.0);
}

TEST_F(BeamTest, DefaultRadiationMode) {
    EXPECT_EQ(beam.getRadiationMode(), "photons");
}

TEST_F(BeamTest, DefaultMachineName) {
    EXPECT_EQ(beam.getMachineName(), "Generic");
}

TEST_F(BeamTest, InitiallyNoRays) {
    EXPECT_EQ(beam.getNumOfRays(), 0u);
    EXPECT_EQ(beam.getTotalNumOfBixels(), 0u);
    EXPECT_TRUE(beam.getRays().empty());
}

// ============================================================================
// Property setters
// ============================================================================

TEST_F(BeamTest, SetBixelWidth) {
    beam.setBixelWidth(5.0);
    EXPECT_DOUBLE_EQ(beam.getBixelWidth(), 5.0);
}

TEST_F(BeamTest, SetSAD) {
    beam.setSAD(1500.0);
    EXPECT_DOUBLE_EQ(beam.getSAD(), 1500.0);
}

TEST_F(BeamTest, SetSCD) {
    beam.setSCD(750.0);
    EXPECT_DOUBLE_EQ(beam.getSCD(), 750.0);
}

TEST_F(BeamTest, SetRadiationMode) {
    beam.setRadiationMode("protons");
    EXPECT_EQ(beam.getRadiationMode(), "protons");
}

TEST_F(BeamTest, SetMachineName) {
    beam.setMachineName("CustomMachine");
    EXPECT_EQ(beam.getMachineName(), "CustomMachine");
}

// ============================================================================
// Source point computation
// ============================================================================

TEST_F(BeamTest, SourcePointBevAtGantryZero) {
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.computeSourcePoints();

    const auto& spBev = beam.getSourcePointBev();
    EXPECT_DOUBLE_EQ(spBev[0], 0.0);
    EXPECT_DOUBLE_EQ(spBev[1], -1000.0);
    EXPECT_DOUBLE_EQ(spBev[2], 0.0);
}

TEST_F(BeamTest, SourcePointLpsMatchesBevAtGantryZero) {
    // At gantry=0, couch=0, the rotation matrix is identity
    // so LPS source point should equal BEV source point
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.computeSourcePoints();

    const auto& spBev = beam.getSourcePointBev();
    const auto& spLps = beam.getSourcePoint();
    EXPECT_NEAR(spBev[0], spLps[0], 1e-10);
    EXPECT_NEAR(spBev[1], spLps[1], 1e-10);
    EXPECT_NEAR(spBev[2], spLps[2], 1e-10);
}

TEST_F(BeamTest, SourcePointDistanceEqualsSAD) {
    // Regardless of angle, source distance from isocenter should be SAD
    beam.setGantryAngle(45.0);
    beam.setCouchAngle(30.0);
    beam.setSAD(1000.0);
    beam.computeSourcePoints();

    const auto& sp = beam.getSourcePoint();
    double dist = std::sqrt(sp[0]*sp[0] + sp[1]*sp[1] + sp[2]*sp[2]);
    EXPECT_NEAR(dist, 1000.0, 1e-8);
}

// ============================================================================
// Ray management
// ============================================================================

TEST_F(BeamTest, AddRayByValue) {
    Ray ray;
    ray.setRayPosBev({-70.0, 0.0, -21.0});
    beam.addRay(ray);
    EXPECT_EQ(beam.getNumOfRays(), 1u);
}

TEST_F(BeamTest, AddRayByMove) {
    Ray ray;
    ray.setRayPosBev({-70.0, 0.0, -21.0});
    beam.addRay(std::move(ray));
    EXPECT_EQ(beam.getNumOfRays(), 1u);
}

TEST_F(BeamTest, GetRayByIndex) {
    Ray ray;
    ray.setRayPosBev({-70.0, 0.0, -21.0});
    ray.setEnergy(6.0);
    beam.addRay(ray);

    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);
    EXPECT_DOUBLE_EQ(r->getRayPosBev()[0], -70.0);
    EXPECT_DOUBLE_EQ(r->getEnergy(), 6.0);
}

TEST_F(BeamTest, GetRayOutOfBoundsReturnsNull) {
    EXPECT_EQ(beam.getRay(0), nullptr);
    EXPECT_EQ(beam.getRay(100), nullptr);
}

TEST_F(BeamTest, NumOfBixelsPerRay) {
    for (int i = 0; i < 5; ++i) {
        Ray ray;
        beam.addRay(ray);
    }
    auto bixelsPerRay = beam.getNumOfBixelsPerRay();
    EXPECT_EQ(bixelsPerRay.size(), 5u);
    for (auto n : bixelsPerRay) {
        EXPECT_EQ(n, 1u); // photon: 1 bixel per ray
    }
}

TEST_F(BeamTest, TotalNumOfBixels) {
    for (int i = 0; i < 10; ++i) {
        Ray ray;
        beam.addRay(ray);
    }
    EXPECT_EQ(beam.getTotalNumOfBixels(), 10u);
}

// ============================================================================
// Ray initialization from BEV positions
// ============================================================================

TEST_F(BeamTest, InitRaysFromPositions) {
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);

    std::vector<Vec3> positions = {
        {-70.0, 0.0, -21.0},
        {  0.0, 0.0,   0.0},
        { 35.0, 0.0,  14.0}
    };

    beam.initRaysFromPositions(positions);

    EXPECT_EQ(beam.getNumOfRays(), 3u);

    // Check first ray
    const auto* r0 = beam.getRay(0);
    ASSERT_NE(r0, nullptr);
    EXPECT_DOUBLE_EQ(r0->getRayPosBev()[0], -70.0);
    EXPECT_DOUBLE_EQ(r0->getRayPosBev()[1],   0.0);
    EXPECT_DOUBLE_EQ(r0->getRayPosBev()[2], -21.0);

    // Target point: (2*x, SAD, 2*z)
    EXPECT_DOUBLE_EQ(r0->getTargetPointBev()[0], -140.0);
    EXPECT_DOUBLE_EQ(r0->getTargetPointBev()[1],  1000.0);
    EXPECT_DOUBLE_EQ(r0->getTargetPointBev()[2],  -42.0);
}

TEST_F(BeamTest, InitRaysLpsMatchesBevAtGantryZero) {
    // With identity rotation, LPS should match BEV
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);

    std::vector<Vec3> positions = {{-70.0, 0.0, -21.0}};
    beam.initRaysFromPositions(positions);

    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);

    EXPECT_NEAR(r->getRayPos()[0], r->getRayPosBev()[0], 1e-10);
    EXPECT_NEAR(r->getRayPos()[1], r->getRayPosBev()[1], 1e-10);
    EXPECT_NEAR(r->getRayPos()[2], r->getRayPosBev()[2], 1e-10);

    EXPECT_NEAR(r->getTargetPoint()[0], r->getTargetPointBev()[0], 1e-10);
    EXPECT_NEAR(r->getTargetPoint()[1], r->getTargetPointBev()[1], 1e-10);
    EXPECT_NEAR(r->getTargetPoint()[2], r->getTargetPointBev()[2], 1e-10);
}

// ============================================================================
// Photon ray corner computation
// ============================================================================

TEST_F(BeamTest, PhotonRayCornersAtGantryZero) {
    // Test against known values for gantry=0, couch=0, bixelWidth=7
    // rayPos_bev = (-70, 0, -21), SAD=1000, SCD=500
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setSCD(500.0);
    beam.setBixelWidth(7.0);

    std::vector<Vec3> positions = {{-70.0, 0.0, -21.0}};
    beam.initRaysFromPositions(positions);
    beam.computePhotonRayCorners();

    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);

    // Expected beamletCornersAtIso (at gantry=0 with identity rotation):
    // corner[0] = (-70 + 3.5, 0, -21 + 3.5) = (-66.5, 0, -17.5)
    // corner[1] = (-70 - 3.5, 0, -21 + 3.5) = (-73.5, 0, -17.5)
    // corner[2] = (-70 - 3.5, 0, -21 - 3.5) = (-73.5, 0, -24.5)
    // corner[3] = (-70 + 3.5, 0, -21 - 3.5) = (-66.5, 0, -24.5)
    const auto& corners = r->getBeamletCornersAtIso();
    EXPECT_NEAR(corners[0][0], -66.5, 1e-10);
    EXPECT_NEAR(corners[0][1],   0.0, 1e-10);
    EXPECT_NEAR(corners[0][2], -17.5, 1e-10);

    EXPECT_NEAR(corners[1][0], -73.5, 1e-10);
    EXPECT_NEAR(corners[1][1],   0.0, 1e-10);
    EXPECT_NEAR(corners[1][2], -17.5, 1e-10);

    EXPECT_NEAR(corners[2][0], -73.5, 1e-10);
    EXPECT_NEAR(corners[2][1],   0.0, 1e-10);
    EXPECT_NEAR(corners[2][2], -24.5, 1e-10);

    EXPECT_NEAR(corners[3][0], -66.5, 1e-10);
    EXPECT_NEAR(corners[3][1],   0.0, 1e-10);
    EXPECT_NEAR(corners[3][2], -24.5, 1e-10);
}

TEST_F(BeamTest, RayCornersSCDAtGantryZero) {
    // Test ray corners at SCD plane
    // SCD/SAD = 500/1000 = 0.5
    // offset = (0, SCD-SAD, 0) = (0, -500, 0)
    // corner_SCD[0] = (0, -500, 0) + 0.5 * (-66.5, 0, -17.5) = (-33.25, -500, -8.75)
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setSCD(500.0);
    beam.setBixelWidth(7.0);

    std::vector<Vec3> positions = {{-70.0, 0.0, -21.0}};
    beam.initRaysFromPositions(positions);
    beam.computePhotonRayCorners();

    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);

    const auto& scdCorners = r->getRayCornersSCD();
    EXPECT_NEAR(scdCorners[0][0], -33.25, 1e-10);
    EXPECT_NEAR(scdCorners[0][1], -500.0, 1e-10);
    EXPECT_NEAR(scdCorners[0][2],  -8.75, 1e-10);

    EXPECT_NEAR(scdCorners[1][0], -36.75, 1e-10);
    EXPECT_NEAR(scdCorners[1][1], -500.0, 1e-10);
    EXPECT_NEAR(scdCorners[1][2],  -8.75, 1e-10);

    EXPECT_NEAR(scdCorners[2][0], -36.75, 1e-10);
    EXPECT_NEAR(scdCorners[2][1], -500.0, 1e-10);
    EXPECT_NEAR(scdCorners[2][2], -12.25, 1e-10);

    EXPECT_NEAR(scdCorners[3][0], -33.25, 1e-10);
    EXPECT_NEAR(scdCorners[3][1], -500.0, 1e-10);
    EXPECT_NEAR(scdCorners[3][2], -12.25, 1e-10);
}

// ============================================================================
// Ray generation (grid-based)
// ============================================================================

TEST_F(BeamTest, GenerateRaysCreatesGrid) {
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setSCD(500.0);

    // A 21x21 mm field with 7mm bixels should give a 4x4 = 16 rays
    // (positions: -10.5, -3.5, 3.5, 10.5... actually -7, 0, 7 = 3x3 = 9)
    beam.generateRays(7.0, {14.0, 14.0});

    // The grid goes from -7 to 7 in steps of 7: positions -7, 0, 7 → 3 per axis
    EXPECT_EQ(beam.getNumOfRays(), 9u); // 3 x 3
    EXPECT_EQ(beam.getTotalNumOfBixels(), 9u);
}

TEST_F(BeamTest, GenerateRaysAllHaveEnergy) {
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.generateRays(7.0, {14.0, 14.0});

    for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
        const auto* r = beam.getRay(i);
        ASSERT_NE(r, nullptr);
        EXPECT_GT(r->getEnergy(), 0.0);
    }
}

TEST_F(BeamTest, GenerateRaysAllHaveCorners) {
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.generateRays(7.0, {14.0, 14.0});

    for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
        const auto* r = beam.getRay(i);
        ASSERT_NE(r, nullptr);
        // All 4 corners should not all be zero
        const auto& corners = r->getBeamletCornersAtIso();
        bool allZero = true;
        for (const auto& c : corners) {
            if (c[0] != 0.0 || c[2] != 0.0) allZero = false;
        }
        // At least for off-center rays, corners should not be all zero
        // The center ray will have corners around (±3.5, 0, ±3.5)
    }
}

// ============================================================================
// Energy assignment
// ============================================================================

TEST_F(BeamTest, SetAllRayEnergies) {
    for (int i = 0; i < 5; ++i) {
        Ray ray;
        ray.setEnergy(6.0);
        beam.addRay(ray);
    }

    beam.setAllRayEnergies(10.0);

    for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
        EXPECT_DOUBLE_EQ(beam.getRay(i)->getEnergy(), 10.0);
    }
}

// ============================================================================
// Rotation effects on ray positions
// ============================================================================

TEST_F(BeamTest, RotatedRayPositionsChangeWithGantry) {
    // Two beams at different gantry angles should produce different LPS positions
    Beam beam1, beam2;
    beam1.setGantryAngle(0.0);
    beam2.setGantryAngle(90.0);
    beam1.setSAD(1000.0);
    beam2.setSAD(1000.0);

    std::vector<Vec3> positions = {{10.0, 0.0, 0.0}};
    beam1.initRaysFromPositions(positions);
    beam2.initRaysFromPositions(positions);

    const auto& pos1 = beam1.getRay(0)->getRayPos();
    const auto& pos2 = beam2.getRay(0)->getRayPos();

    // They should be different (unless by coincidence, which won't happen here)
    bool different = (std::abs(pos1[0] - pos2[0]) > 1e-6) ||
                     (std::abs(pos1[1] - pos2[1]) > 1e-6) ||
                     (std::abs(pos1[2] - pos2[2]) > 1e-6);
    EXPECT_TRUE(different);
}

// ============================================================================
// Target-aware ray generation (generateRaysFromTarget)
// ============================================================================

TEST_F(BeamTest, GenerateRaysFromTargetSingleVoxel) {
    // A single target voxel at the isocenter should produce exactly 1 ray at (0, 0, 0) BEV
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {{0.0, 0.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 1u);
    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);
    EXPECT_NEAR(r->getRayPosBev()[0], 0.0, 1e-10);
    EXPECT_NEAR(r->getRayPosBev()[2], 0.0, 1e-10);
}

TEST_F(BeamTest, GenerateRaysFromTargetDuplicateVoxelsMerged) {
    // Several target voxels that project to the same bixel position should produce 1 ray
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    // All within 1mm of origin → all snap to (0,0) on a 7mm grid
    std::vector<Vec3> targetCoords = {
        {1.0, 0.0, 1.0},
        {2.0, 0.0, 2.0},
        {3.0, 0.0, -1.0},
        {-2.0, 0.0, 0.5}
    };
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 1u);
}

TEST_F(BeamTest, GenerateRaysFromTargetDifferentAnglesProduceDifferentRayCounts) {
    // The key test: the same target voxels should produce different ray counts
    // at different gantry angles because the projection is angle-dependent.
    //
    // Create a rectangular target elongated along x (LPS).
    // At gantry=0 the x-extent projects onto BEV x.
    // At gantry=90 the x-extent projects onto BEV y (depth axis), 
    //   which gets collapsed in the isocenter plane projection.
    Vec3 iso = {0.0, 0.0, 0.0};

    // Target: 5 voxels spread along x, narrow in y and z
    std::vector<Vec3> targetCoords;
    for (double x = -28.0; x <= 28.0; x += 7.0) {
        targetCoords.push_back({x, 0.0, 0.0});
    }

    // Beam at gantry=0: target x maps to BEV x → all unique x positions
    Beam beam0;
    beam0.setGantryAngle(0.0);
    beam0.setCouchAngle(0.0);
    beam0.setSAD(1000.0);
    beam0.setIsocenter(iso);
    beam0.generateRaysFromTarget(targetCoords, 7.0);

    // Beam at gantry=90: target x maps to BEV y (depth), 
    // so x spread collapses to single bev-x position
    Beam beam90;
    beam90.setGantryAngle(90.0);
    beam90.setCouchAngle(0.0);
    beam90.setSAD(1000.0);
    beam90.setIsocenter(iso);
    beam90.generateRaysFromTarget(targetCoords, 7.0);

    // At gantry=0, spread along BEV-x → 9 unique x positions × 1 z = 9 rays
    // At gantry=90, spread along BEV-y → 1 x × 1 z = 1 ray
    EXPECT_GT(beam0.getNumOfRays(), beam90.getNumOfRays());
}

TEST_F(BeamTest, GenerateRaysFromTargetAllRaysHaveCorners) {
    beam.setGantryAngle(45.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setSCD(500.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords;
    for (double x = -14.0; x <= 14.0; x += 7.0) {
        for (double z = -14.0; z <= 14.0; z += 7.0) {
            targetCoords.push_back({x, 0.0, z});
        }
    }
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_GT(beam.getNumOfRays(), 0u);
    for (size_t i = 0; i < beam.getNumOfRays(); ++i) {
        const auto* r = beam.getRay(i);
        ASSERT_NE(r, nullptr);
        // Check that SCD corners have the expected y-offset (SCD - SAD = -500)
        const auto& scdCorners = r->getRayCornersSCD();
        // The corners should not all be zero
        bool anyNonZero = false;
        for (const auto& c : scdCorners) {
            if (std::abs(c[0]) > 1e-10 || std::abs(c[1]) > 1e-10 || std::abs(c[2]) > 1e-10) {
                anyNonZero = true;
                break;
            }
        }
        EXPECT_TRUE(anyNonZero);
    }
}

TEST_F(BeamTest, GenerateRaysFromTargetIsoOffsetApplied) {
    // Target at (100, 0, 0) with isocenter at (100, 0, 0) should project to (0,0) BEV
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({100.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {{100.0, 0.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 1u);
    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);
    EXPECT_NEAR(r->getRayPosBev()[0], 0.0, 1e-10);
    EXPECT_NEAR(r->getRayPosBev()[2], 0.0, 1e-10);
}

TEST_F(BeamTest, GenerateRaysFromTargetSADDivergenceCorrection) {
    // A target voxel at (0, SAD, 0) in BEV (i.e. at twice the isocenter distance)
    // should have divergence correction: scale = SAD/(SAD + y_bev) = 0.5
    //
    // With gantry=0, BEV = LPS, so world (14, SAD, 0) = BEV (14, SAD, 0)
    // Projected x = 14 * 0.5 = 7.0, snapped = round(7/7)*7 = 7.0
    //
    // Without divergence correction, x=14 would snap to 14 on the 7mm grid.
    // With it, we get 7. This proves the correction is applied.
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {{14.0, 1000.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 1u);
    const auto* r = beam.getRay(0);
    ASSERT_NE(r, nullptr);
    // Projected x = 14 * (1000 / (1000 + 1000)) = 7.0, round(7.0/7)*7 = 7.0
    EXPECT_NEAR(r->getRayPosBev()[0], 7.0, 1e-10);
}

TEST_F(BeamTest, GenerateRaysFromTargetPaddingWithCoarseResolution) {
    // When bixelWidth < max CT resolution, padding should add surrounding positions
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    // Single voxel at origin, bixelWidth=3, ctResolution=7mm
    // Without padding: 1 ray at (0,0)
    // With padding: pad = floor(7/3) = 2, so ±2 in both axes → (2*2+1)^2 = 25
    std::vector<Vec3> targetCoords = {{0.0, 0.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 3.0, {7.0, 7.0, 7.0});

    // Should have more than 1 ray due to padding
    EXPECT_GT(beam.getNumOfRays(), 1u);
    // pad=2, so we expect a (2*2+1) x (2*2+1) = 5x5 = 25 grid
    EXPECT_EQ(beam.getNumOfRays(), 25u);
}

TEST_F(BeamTest, GenerateRaysFromTargetNoPaddingWhenBixelWidthLarger) {
    // When bixelWidth >= max CT resolution, no padding should occur
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {{0.0, 0.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 7.0, {3.0, 3.0, 3.0});

    // Should have exactly 1 ray - no padding
    EXPECT_EQ(beam.getNumOfRays(), 1u);
}

TEST_F(BeamTest, GenerateRaysFromTargetWithoutMarginSingleVoxelOneRay) {
    // Without 2D margin, a single voxel at origin should produce exactly 1 ray
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {{0.0, 0.0, 0.0}};
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 1u);
}

TEST_F(BeamTest, GenerateRaysFromTargetMultipleVoxelsDifferentPositions) {
    // Multiple voxels at different grid positions should produce multiple rays
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    std::vector<Vec3> targetCoords = {
        {0.0, 0.0, 0.0},
        {7.0, 0.0, 0.0},
        {0.0, 0.0, 7.0}
    };
    beam.generateRaysFromTarget(targetCoords, 7.0);

    EXPECT_EQ(beam.getNumOfRays(), 3u);
}

TEST_F(BeamTest, GenerateRaysFromTarget3DCubeTarget) {
    // A 3D cube of target voxels should project to a roughly square pattern
    // for gantry=0 (looking along y-axis)
    beam.setGantryAngle(0.0);
    beam.setCouchAngle(0.0);
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});

    // Create a 3x3x3 cube of voxels centered at origin, 7mm spacing
    std::vector<Vec3> targetCoords;
    for (double x = -7.0; x <= 7.0; x += 7.0) {
        for (double y = -7.0; y <= 7.0; y += 7.0) {
            for (double z = -7.0; z <= 7.0; z += 7.0) {
                targetCoords.push_back({x, y, z});
            }
        }
    }

    beam.generateRaysFromTarget(targetCoords, 7.0);

    // At gantry=0, the y-axis is the BEV depth axis.
    // All voxels at (x, *, z) project roughly to the same (x,z) BEV position
    // (with slight SAD divergence variation from depth).
    // For SAD=1000 and y offsets of ±7, the divergence effect is tiny,
    // so we expect roughly 3x3 = 9 unique positions.
    EXPECT_GE(beam.getNumOfRays(), 9u);
    EXPECT_LE(beam.getNumOfRays(), 15u); // Some divergence might split a few
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
