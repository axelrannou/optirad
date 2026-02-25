#include <gtest/gtest.h>
#include "dose/RadDepthCalculator.hpp"
#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include <cmath>

namespace optirad::tests {

class RadDepthCalculatorTest : public ::testing::Test {
protected:
    Grid grid;
    std::vector<double> densityData;

    void SetUp() override {
        // 5x5x5 grid, spacing 10mm, origin at (0,0,0)
        grid.setDimensions(5, 5, 5);
        grid.setSpacing(10.0, 10.0, 10.0);
        grid.setOrigin(Vec3{0.0, 0.0, 0.0});

        // Uniform water density
        densityData.assign(grid.getNumVoxels(), 1.0);
    }
};

// ============================================================================
// computeRadDepths
// ============================================================================

TEST_F(RadDepthCalculatorTest, AllVoxelsHavePositiveDepth) {
    Vec3 source = {-50.0, 25.0, 25.0};

    // Select a few voxels along the x-axis (center row)
    std::vector<size_t> indices;
    auto dims = grid.getDimensions();
    for (size_t ix = 0; ix < 5; ++ix) {
        size_t idx = ix + 2 * dims[0] + 2 * dims[0] * dims[1];
        indices.push_back(idx);
    }

    auto depthMap = RadDepthCalculator::computeRadDepths(
        source, grid, densityData.data(), indices);

    // All voxels should have positive radiological depth (they are in front of source)
    for (size_t ix = 0; ix < 5; ++ix) {
        size_t idx = ix + 2 * dims[0] + 2 * dims[0] * dims[1];
        auto it = depthMap.find(idx);
        ASSERT_NE(it, depthMap.end()) << "Voxel " << idx << " not found in depth map";
        EXPECT_GT(it->second, 0.0) << "Depth should be positive for voxel ix=" << ix;
    }
}

TEST_F(RadDepthCalculatorTest, DoubleDensityDoublesDepth) {
    Vec3 source = {-50.0, 25.0, 25.0};

    // Center voxel
    auto dims = grid.getDimensions();
    size_t centerIdx = 2 + 2 * dims[0] + 2 * dims[0] * dims[1];
    std::vector<size_t> indices = {centerIdx};

    // With density 1.0
    auto depth1 = RadDepthCalculator::computeRadDepths(
        source, grid, densityData.data(), indices);

    // Double density
    std::vector<double> densityX2(densityData.size(), 2.0);
    auto depth2 = RadDepthCalculator::computeRadDepths(
        source, grid, densityX2.data(), indices);

    ASSERT_NE(depth1.find(centerIdx), depth1.end());
    ASSERT_NE(depth2.find(centerIdx), depth2.end());

    // Depth with 2x density should be ~2x depth with 1x density
    EXPECT_NEAR(depth2.at(centerIdx), 2.0 * depth1.at(centerIdx), 
                0.1 * depth1.at(centerIdx));
}

// ============================================================================
// computeRayRadDepths
// ============================================================================

TEST_F(RadDepthCalculatorTest, RayRadDepthsCumulative) {
    Vec3 source = {-50.0, 25.0, 25.0};
    Vec3 target = {75.0, 25.0, 25.0};

    auto depths = RadDepthCalculator::computeRayRadDepths(
        source, target, grid, densityData.data());

    // Should have 5 voxels (along x with y=2, z=2)
    EXPECT_GE(depths.size(), 4u);

    // Cumulative: each depth ≥ previous
    double prev = 0.0;
    for (const auto& [idx, depth] : depths) {
        EXPECT_GE(depth, prev);
        prev = depth;
    }
}

TEST_F(RadDepthCalculatorTest, EmptyVoxelListReturnsEmpty) {
    Vec3 source = {-50.0, 25.0, 25.0};
    std::vector<size_t> indices;

    auto depthMap = RadDepthCalculator::computeRadDepths(
        source, grid, densityData.data(), indices);

    EXPECT_TRUE(depthMap.empty());
}

} // namespace optirad::tests
