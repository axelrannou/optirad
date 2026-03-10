#include <gtest/gtest.h>
#include "dose/SSDCalculator.hpp"
#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include <cmath>

namespace optirad::tests {

class SSDCalculatorTest : public ::testing::Test {
protected:
    Grid grid;
    std::vector<double> densityData;

    void SetUp() override {
        // 10x10x10 grid, spacing 5mm, origin at (0,0,0)
        grid.setDimensions(10, 10, 10);
        grid.setSpacing(5.0, 5.0, 5.0);
        grid.setOrigin(Vec3{0.0, 0.0, 0.0});

        densityData.assign(grid.getNumVoxels(), 0.0);

        // Place a "body" in the center region: voxels [2..7] in each dimension → density = 1.0
        auto dims = grid.getDimensions();
        for (size_t z = 2; z < 8; ++z) {
            for (size_t y = 2; y < 8; ++y) {
                for (size_t x = 2; x < 8; ++x) {
                    size_t idx = x + y * dims[0] + z * dims[0] * dims[1];
                    densityData[idx] = 1.0;
                }
            }
        }
    }
};

// ============================================================================
// Single ray SSD
// ============================================================================

TEST_F(SSDCalculatorTest, RayHittingBody) {
    // Source far left, target in center
    Vec3 source = {-100.0, 25.0, 25.0};
    Vec3 target = {25.0, 25.0, 25.0};

    double ssd = SSDCalculator::computeSSD(source, target, grid, densityData.data());

    // Body surface at x = 10mm (voxel 2 starts at 2*5=10),
    // source at x=-100, so SSD ~110mm
    EXPECT_GT(ssd, 0.0);
    EXPECT_LT(ssd, 200.0);
}

TEST_F(SSDCalculatorTest, RayMissingBody) {
    // Source and target completely outside the body region
    Vec3 source = {-100.0, 0.0, 0.0};
    Vec3 target = {-50.0, 0.0, 0.0};

    double ssd = SSDCalculator::computeSSD(source, target, grid, densityData.data());

    // Ray through air only → SSD = -1 (miss)
    EXPECT_EQ(ssd, -1.0);
}

// ============================================================================
// Beam SSDs
// ============================================================================

TEST_F(SSDCalculatorTest, BeamSSDs) {
    Vec3 source = {-100.0, 25.0, 25.0};

    // Two rays: one hitting center, one at the body corner, both should find surface
    std::vector<Vec3> targets = {
        {25.0, 25.0, 25.0},
        {25.0, 15.0, 25.0}
    };

    auto ssds = SSDCalculator::computeBeamSSDs(
        source, targets, grid, densityData.data());

    ASSERT_EQ(ssds.size(), 2u);
    EXPECT_GT(ssds[0], 0.0);
    EXPECT_GT(ssds[1], 0.0);
}

// ============================================================================
// Density threshold
// ============================================================================

TEST_F(SSDCalculatorTest, HighThresholdMisses) {
    Vec3 source = {-100.0, 25.0, 25.0};
    Vec3 target = {25.0, 25.0, 25.0};

    // With threshold = 2.0, the body (density 1.0) should not count as surface
    double ssd = SSDCalculator::computeSSD(source, target, grid, densityData.data(), 2.0);
    EXPECT_EQ(ssd, -1.0);
}

} // namespace optirad::tests
