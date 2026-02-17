#include <gtest/gtest.h>
#include "geometry/VoxelDilation.hpp"
#include <set>
#include <algorithm>

namespace optirad::tests {

class VoxelDilationTest : public ::testing::Test {
protected:
    // Helper: convert (row, col, slice) to 1-based linear index
    // Column-major: idx1 = row + col*ny + slice*ny*nx + 1
    size_t sub2ind(size_t row, size_t col, size_t slice, size_t ny, size_t nx) const {
        return row + col * ny + slice * ny * nx + 1;
    }

    // Helper: convert 1-based linear index to subscripts
    void ind2sub(size_t idx1, size_t ny, size_t nx, size_t& row, size_t& col, size_t& slice) const {
        size_t idx0 = idx1 - 1;
        row   = idx0 % ny;
        col   = (idx0 / ny) % nx;
        slice = idx0 / (ny * nx);
    }
};

TEST_F(VoxelDilationTest, NoMarginReturnsOriginal) {
    std::array<size_t, 3> dims = {10, 10, 10};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {0.0, 0.0, 0.0};

    std::vector<size_t> target = {sub2ind(5, 5, 5, 10, 10)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 1000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);
    EXPECT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0], target[0]);
}

TEST_F(VoxelDilationTest, SingleVoxelDilatedBy1) {
    // 10x10x10 grid, single voxel at center (5,5,5), margin=1mm, spacing=1mm
    // Should expand to 3x3x3 = 27 voxels (26-connectivity, 1 iteration)
    std::array<size_t, 3> dims = {10, 10, 10};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {1.0, 1.0, 1.0};

    std::vector<size_t> target = {sub2ind(5, 5, 5, 10, 10)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 1000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // 26-connectivity dilation of a single voxel by 1 step: 3x3x3 = 27
    EXPECT_EQ(result.size(), 27u);
}

TEST_F(VoxelDilationTest, ConstrainedToPatientSurface) {
    // Only include a subset of voxels as "patient surface"
    // Dilation should not extend beyond patient surface
    std::array<size_t, 3> dims = {10, 10, 10};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {1.0, 1.0, 1.0};

    std::vector<size_t> target = {sub2ind(5, 5, 5, 10, 10)};

    // Only add a few neighbors as patient surface
    std::set<size_t> allVoxels;
    allVoxels.insert(sub2ind(5, 5, 5, 10, 10)); // center
    allVoxels.insert(sub2ind(6, 5, 5, 10, 10)); // +row
    allVoxels.insert(sub2ind(4, 5, 5, 10, 10)); // -row
    allVoxels.insert(sub2ind(5, 6, 5, 10, 10)); // +col

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // Should only include voxels that are in allVoxels
    std::set<size_t> resultSet(result.begin(), result.end());
    for (size_t idx : resultSet) {
        EXPECT_TRUE(allVoxels.count(idx) > 0);
    }
    // Should have at most 4 voxels
    EXPECT_LE(result.size(), 4u);
}

TEST_F(VoxelDilationTest, BorderVoxelsNotExpanded) {
    // Voxel at the border of the grid should not cause out-of-bounds
    std::array<size_t, 3> dims = {10, 10, 10};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {1.0, 1.0, 1.0};

    // Voxel at corner (0, 0, 0)
    std::vector<size_t> target = {sub2ind(0, 0, 0, 10, 10)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 1000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // Border voxels are excluded from expansion, so only the original voxel
    EXPECT_GE(result.size(), 1u);
}

TEST_F(VoxelDilationTest, AnisotropicMargin) {
    // Different margin per axis: margin_y=2, margin_x=1, margin_z=0
    // voxelMargins = round({2/1, 1/1, 0/1}) = {2, 1, 0}
    // maxIter = 2
    std::array<size_t, 3> dims = {20, 20, 20};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {2.0, 1.0, 0.0};

    std::vector<size_t> target = {sub2ind(10, 10, 10, 20, 20)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 8000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // The expansion should be asymmetric
    // Row (y): ±2, Col (x): ±1, Slice (z): 0
    std::set<size_t> resultSet(result.begin(), result.end());

    // Check center is included
    EXPECT_TRUE(resultSet.count(sub2ind(10, 10, 10, 20, 20)));

    // Check row expansion: (12, 10, 10) should be included (2 steps in row)
    EXPECT_TRUE(resultSet.count(sub2ind(12, 10, 10, 20, 20)));

    // Check col expansion: (10, 11, 10) should be included (1 step in col)
    EXPECT_TRUE(resultSet.count(sub2ind(10, 11, 10, 20, 20)));

    // Check no slice expansion: (10, 10, 11) should NOT be included
    EXPECT_FALSE(resultSet.count(sub2ind(10, 10, 11, 20, 20)));
}

TEST_F(VoxelDilationTest, MultipleIterationsExpand) {
    // Margin = 2mm, spacing = 1mm → voxelMargins = {2,2,2}, maxIter = 2
    // Single voxel should expand to 5x5x5 = 125
    std::array<size_t, 3> dims = {20, 20, 20};
    std::array<double, 3> spacing = {1.0, 1.0, 1.0};
    std::array<double, 3> margin = {2.0, 2.0, 2.0};

    std::vector<size_t> target = {sub2ind(10, 10, 10, 20, 20)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 8000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // 26-connectivity dilation by 2: in each iteration, each active voxel
    // expands to its 26 neighbors. After 2 iterations from a single voxel,
    // the result is a 5x5x5 cube = 125 voxels
    EXPECT_EQ(result.size(), 125u);
}

TEST_F(VoxelDilationTest, NonUniformSpacingCorrectVoxelMargins) {
    // margin = {7, 7, 7} mm, spacing = {0.976, 0.976, 2.0} mm
    // voxelMargins = round({7/0.976, 7/0.976, 7/2.0}) = round({7.17, 7.17, 3.5}) = {7, 7, 4}
    // This mimics the JOHN_DOE_LUNG case
    std::array<size_t, 3> dims = {20, 20, 20};
    std::array<double, 3> spacing = {0.976562, 0.976562, 2.0};
    std::array<double, 3> margin = {7.0, 7.0, 7.0};

    std::vector<size_t> target = {sub2ind(10, 10, 10, 20, 20)};
    std::set<size_t> allVoxels;
    for (size_t i = 1; i <= 8000; ++i) allVoxels.insert(i);

    auto result = dilateVoxels(target, allVoxels, dims, spacing, margin);

    // voxelMargins = {7, 7, 4}, maxIter = 7
    // Result should be (2*7+1) × (2*7+1) × (2*4+1) = 15 × 15 × 9 = 2025 voxels
    // but limited by grid boundaries (center at 10, grid 0-19)
    // In row: 10±7 = 3..17 → 15 voxels (fits)
    // In col: 10±7 = 3..17 → 15 voxels (fits)
    // In slice: 10±4 = 6..14 → 9 voxels (fits)
    // So 15 * 15 * 9 = 2025
    EXPECT_EQ(result.size(), 2025u);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
