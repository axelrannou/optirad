#include <gtest/gtest.h>
#include "core/Ray.hpp"
#include <cmath>

namespace optirad::tests {

class RayTest : public ::testing::Test {
protected:
    Ray ray;
};

// --- Default construction ---

TEST_F(RayTest, DefaultRayPosBevIsZero) {
    const auto& pos = ray.getRayPosBev();
    EXPECT_DOUBLE_EQ(pos[0], 0.0);
    EXPECT_DOUBLE_EQ(pos[1], 0.0);
    EXPECT_DOUBLE_EQ(pos[2], 0.0);
}

TEST_F(RayTest, DefaultTargetPointBevIsZero) {
    const auto& tp = ray.getTargetPointBev();
    EXPECT_DOUBLE_EQ(tp[0], 0.0);
    EXPECT_DOUBLE_EQ(tp[1], 0.0);
    EXPECT_DOUBLE_EQ(tp[2], 0.0);
}

TEST_F(RayTest, DefaultRayPosIsZero) {
    const auto& pos = ray.getRayPos();
    EXPECT_DOUBLE_EQ(pos[0], 0.0);
    EXPECT_DOUBLE_EQ(pos[1], 0.0);
    EXPECT_DOUBLE_EQ(pos[2], 0.0);
}

TEST_F(RayTest, DefaultEnergyIs6MV) {
    EXPECT_DOUBLE_EQ(ray.getEnergy(), 6.0);
}

TEST_F(RayTest, DefaultNumOfBixelsIs1) {
    EXPECT_EQ(ray.getNumOfBixels(), 1u);
}

// --- Setters and getters ---

TEST_F(RayTest, SetRayPosBev) {
    ray.setRayPosBev({-70.0, 0.0, -21.0});
    const auto& pos = ray.getRayPosBev();
    EXPECT_DOUBLE_EQ(pos[0], -70.0);
    EXPECT_DOUBLE_EQ(pos[1], 0.0);
    EXPECT_DOUBLE_EQ(pos[2], -21.0);
}

TEST_F(RayTest, SetTargetPointBev) {
    ray.setTargetPointBev({-140.0, 1000.0, -42.0});
    const auto& tp = ray.getTargetPointBev();
    EXPECT_DOUBLE_EQ(tp[0], -140.0);
    EXPECT_DOUBLE_EQ(tp[1], 1000.0);
    EXPECT_DOUBLE_EQ(tp[2], -42.0);
}

TEST_F(RayTest, SetRayPos) {
    ray.setRayPos({10.0, 20.0, 30.0});
    const auto& pos = ray.getRayPos();
    EXPECT_DOUBLE_EQ(pos[0], 10.0);
    EXPECT_DOUBLE_EQ(pos[1], 20.0);
    EXPECT_DOUBLE_EQ(pos[2], 30.0);
}

TEST_F(RayTest, SetTargetPoint) {
    ray.setTargetPoint({100.0, 200.0, 300.0});
    const auto& tp = ray.getTargetPoint();
    EXPECT_DOUBLE_EQ(tp[0], 100.0);
    EXPECT_DOUBLE_EQ(tp[1], 200.0);
    EXPECT_DOUBLE_EQ(tp[2], 300.0);
}

TEST_F(RayTest, SetEnergy) {
    ray.setEnergy(10.0);
    EXPECT_DOUBLE_EQ(ray.getEnergy(), 10.0);
}

// --- Beamlet corners ---

TEST_F(RayTest, SetBeamletCornersAtIso) {
    std::array<Vec3, 4> corners = {{
        {-66.5, 0.0, -17.5},
        {-73.5, 0.0, -17.5},
        {-73.5, 0.0, -24.5},
        {-66.5, 0.0, -24.5}
    }};
    ray.setBeamletCornersAtIso(corners);

    const auto& result = ray.getBeamletCornersAtIso();
    EXPECT_DOUBLE_EQ(result[0][0], -66.5);
    EXPECT_DOUBLE_EQ(result[1][0], -73.5);
    EXPECT_DOUBLE_EQ(result[2][2], -24.5);
    EXPECT_DOUBLE_EQ(result[3][2], -24.5);
}

TEST_F(RayTest, SetRayCornersSCD) {
    std::array<Vec3, 4> corners = {{
        {-33.25, -500.0, -8.75},
        {-36.75, -500.0, -8.75},
        {-36.75, -500.0, -12.25},
        {-33.25, -500.0, -12.25}
    }};
    ray.setRayCornersSCD(corners);

    const auto& result = ray.getRayCornersSCD();
    EXPECT_DOUBLE_EQ(result[0][0], -33.25);
    EXPECT_DOUBLE_EQ(result[0][1], -500.0);
    EXPECT_DOUBLE_EQ(result[3][2], -12.25);
}

// --- Consistency between BEV and LPS for zero gantry/couch ---

TEST_F(RayTest, BevAndLpsConsistentForIdentityRotation) {
    // For gantry=0, couch=0, BEV and LPS positions should match
    // (since rotation matrix is identity)
    Vec3 pos = {-70.0, 0.0, -21.0};
    Vec3 target = {-140.0, 1000.0, -42.0};
    
    ray.setRayPosBev(pos);
    ray.setRayPos(pos);  // Would be set by Beam after rotation
    ray.setTargetPointBev(target);
    ray.setTargetPoint(target);

    EXPECT_DOUBLE_EQ(ray.getRayPosBev()[0], ray.getRayPos()[0]);
    EXPECT_DOUBLE_EQ(ray.getRayPosBev()[1], ray.getRayPos()[1]);
    EXPECT_DOUBLE_EQ(ray.getRayPosBev()[2], ray.getRayPos()[2]);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
