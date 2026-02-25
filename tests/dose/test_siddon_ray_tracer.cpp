#include <gtest/gtest.h>
#include "dose/SiddonRayTracer.hpp"
#include "geometry/Grid.hpp"
#include "geometry/MathUtils.hpp"
#include <cmath>

namespace optirad::tests {

class SiddonRayTracerTest : public ::testing::Test {
protected:
    Grid grid;
    std::vector<double> densityData;

    void SetUp() override {
        // Simple 5x5x5 grid, origin at (0,0,0), spacing 10mm
        grid.setDimensions(5, 5, 5);
        grid.setSpacing(10.0, 10.0, 10.0);
        grid.setOrigin(Vec3{0.0, 0.0, 0.0});

        // Uniform water density (1.0)
        densityData.assign(grid.getNumVoxels(), 1.0);
    }
};

// ============================================================================
// Ray along principal axis
// ============================================================================

TEST_F(SiddonRayTracerTest, RayAlongXAxis) {
    // Ray from left of grid to right, through center row/col
    Vec3 source = {-20.0, 25.0, 25.0}; // Outside left
    Vec3 target = {70.0, 25.0, 25.0};  // Outside right

    auto result = SiddonRayTracer::trace(source, target, grid, densityData.data());

    // Should traverse 5 voxels along x
    EXPECT_EQ(result.voxelIndices.size(), 5u);
    EXPECT_GT(result.totalDistance, 0.0);

    // Each voxel should have ~10mm intersection length
    for (double len : result.intersectionLengths) {
        EXPECT_NEAR(len, 10.0, 0.5);
    }
}

TEST_F(SiddonRayTracerTest, RayAlongZAxis) {
    Vec3 source = {25.0, 25.0, -20.0};
    Vec3 target = {25.0, 25.0, 70.0};

    auto result = SiddonRayTracer::trace(source, target, grid, densityData.data());

    EXPECT_EQ(result.voxelIndices.size(), 5u);
    for (double len : result.intersectionLengths) {
        EXPECT_NEAR(len, 10.0, 0.5);
    }
}

// ============================================================================
// Ray misses the volume
// ============================================================================

TEST_F(SiddonRayTracerTest, RayMissesVolume) {
    Vec3 source = {-20.0, 100.0, 100.0}; // Far from grid
    Vec3 target = {70.0, 100.0, 100.0};

    auto result = SiddonRayTracer::trace(source, target, grid, densityData.data());

    EXPECT_TRUE(result.voxelIndices.empty());
}

// ============================================================================
// Diagonal ray traverses more voxels
// ============================================================================

TEST_F(SiddonRayTracerTest, DiagonalRayTraversesMultipleVoxels) {
    Vec3 source = {-10.0, -10.0, 25.0};
    Vec3 target = {60.0, 60.0, 25.0};

    auto result = SiddonRayTracer::trace(source, target, grid, densityData.data());

    // Diagonal: expect traversal through multiple voxels (at least 5)
    EXPECT_GE(result.voxelIndices.size(), 5u);
}

// ============================================================================
// traceRadDepth
// ============================================================================

TEST_F(SiddonRayTracerTest, RadDepthAccumulates) {
    Vec3 source = {-20.0, 25.0, 25.0};
    Vec3 target = {70.0, 25.0, 25.0};

    auto radDepths = SiddonRayTracer::traceRadDepth(
        source, target, grid, densityData.data());

    // Should have entries and be cumulative (increasing)
    EXPECT_FALSE(radDepths.empty());

    double prevDepth = 0.0;
    for (const auto& [idx, depth] : radDepths) {
        EXPECT_GE(depth, prevDepth);
        prevDepth = depth;
    }

    // With unit density and 10mm spacing, total radDepth ≈ 50mm
    if (!radDepths.empty()) {
        EXPECT_NEAR(radDepths.back().second, 50.0, 5.0);
    }
}

// ============================================================================
// Null density data
// ============================================================================

TEST_F(SiddonRayTracerTest, TraceWithNullDensity) {
    Vec3 source = {-20.0, 25.0, 25.0};
    Vec3 target = {70.0, 25.0, 25.0};

    auto result = SiddonRayTracer::trace(source, target, grid, nullptr);

    // Should still trace geometry
    EXPECT_GT(result.voxelIndices.size(), 0u);
}

// ============================================================================
// Variable density
// ============================================================================

TEST_F(SiddonRayTracerTest, RadDepthScalesWithDensity) {
    // Set density to 2.0 everywhere
    std::fill(densityData.begin(), densityData.end(), 2.0);

    Vec3 source = {-20.0, 25.0, 25.0};
    Vec3 target = {70.0, 25.0, 25.0};

    auto radDepths = SiddonRayTracer::traceRadDepth(
        source, target, grid, densityData.data());

    // With density 2.0 everywhere, rad depth should be ~100mm (50mm * 2.0)
    if (!radDepths.empty()) {
        EXPECT_NEAR(radDepths.back().second, 100.0, 10.0);
    }
}

} // namespace optirad::tests
