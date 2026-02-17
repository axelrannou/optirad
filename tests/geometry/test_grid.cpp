#include <gtest/gtest.h>
#include "geometry/Grid.hpp"

namespace optirad::tests {

class GridTest : public ::testing::Test {
protected:
    Grid grid;
    
    void SetUp() override {
        // Create a 10x10x10 grid with 1mm spacing
        grid.setDimensions(10, 10, 10);
        grid.setSpacing(1.0, 1.0, 1.0);
        grid.setOrigin({0.0, 0.0, 0.0});
    }
};

TEST_F(GridTest, SetDimensionsWorks) {
    auto dims = grid.getDimensions();
    EXPECT_EQ(dims[0], 10);
    EXPECT_EQ(dims[1], 10);
    EXPECT_EQ(dims[2], 10);
}

TEST_F(GridTest, SetSpacingWorks) {
    auto spacing = grid.getSpacing();
    EXPECT_DOUBLE_EQ(spacing[0], 1.0);
    EXPECT_DOUBLE_EQ(spacing[1], 1.0);
    EXPECT_DOUBLE_EQ(spacing[2], 1.0);
}

TEST_F(GridTest, SetOriginWorks) {
    auto origin = grid.getOrigin();
    EXPECT_DOUBLE_EQ(origin[0], 0.0);
    EXPECT_DOUBLE_EQ(origin[1], 0.0);
    EXPECT_DOUBLE_EQ(origin[2], 0.0);
}

TEST_F(GridTest, GetNumVoxelsIsCorrect) {
    size_t expectedVoxels = 10 * 10 * 10;
    EXPECT_EQ(grid.getNumVoxels(), expectedVoxels);
}

TEST_F(GridTest, NonUniformSpacing) {
    Grid nonUniform;
    nonUniform.setDimensions(5, 5, 5);
    nonUniform.setSpacing(0.5, 1.0, 2.0);
    nonUniform.setOrigin({10.0, 20.0, 30.0});
    
    auto spacing = nonUniform.getSpacing();
    EXPECT_DOUBLE_EQ(spacing[0], 0.5);
    EXPECT_DOUBLE_EQ(spacing[1], 1.0);
    EXPECT_DOUBLE_EQ(spacing[2], 2.0);
    
    auto origin = nonUniform.getOrigin();
    EXPECT_DOUBLE_EQ(origin[0], 10.0);
    EXPECT_DOUBLE_EQ(origin[1], 20.0);
    EXPECT_DOUBLE_EQ(origin[2], 30.0);
}

TEST_F(GridTest, LargeGridDimensions) {
    Grid large;
    large.setDimensions(512, 512, 256);
    auto dims = large.getDimensions();
    EXPECT_EQ(dims[0], 512);
    EXPECT_EQ(dims[1], 512);
    EXPECT_EQ(dims[2], 256);
}

TEST_F(GridTest, PatientPositionProperty) {
    grid.setPatientPosition("HFS");
    EXPECT_EQ(grid.getPatientPosition(), "HFS");
}

TEST_F(GridTest, SliceThicknessProperty) {
    grid.setSliceThickness(5.0);
    EXPECT_DOUBLE_EQ(grid.getSliceThickness(), 5.0);
}

TEST_F(GridTest, CoordinateArraysHaveCorrectSize) {
    auto x = grid.getXCoordinates();
    auto y = grid.getYCoordinates();
    auto z = grid.getZCoordinates();
    
    EXPECT_EQ(x.size(), 10);
    EXPECT_EQ(y.size(), 10);
    EXPECT_EQ(z.size(), 10);
}

TEST_F(GridTest, CoordinateArraysWithOriginZero) {
    // With origin at (0,0,0) and unit spacing
    auto x = grid.getXCoordinates();
    auto y = grid.getYCoordinates();
    auto z = grid.getZCoordinates();
    
    // First voxel center should be at origin
    EXPECT_DOUBLE_EQ(x[0], 0.0);
    EXPECT_DOUBLE_EQ(y[0], 0.0);
    EXPECT_DOUBLE_EQ(z[0], 0.0);
    
    // Last voxel center should be at (dimensions-1) * spacing
    EXPECT_DOUBLE_EQ(x[9], 9.0);
    EXPECT_DOUBLE_EQ(y[9], 9.0);
    EXPECT_DOUBLE_EQ(z[9], 9.0);
}

TEST_F(GridTest, CoordinateArraysWithNonZeroOrigin) {
    Grid shifted;
    shifted.setDimensions(5, 5, 5);
    shifted.setSpacing(2.0, 2.0, 2.0);
    shifted.setOrigin({-10.0, -20.0, -30.0});
    
    auto x = shifted.getXCoordinates();
    auto y = shifted.getYCoordinates();
    auto z = shifted.getZCoordinates();
    
    // First voxel should be at origin
    EXPECT_DOUBLE_EQ(x[0], -10.0);
    EXPECT_DOUBLE_EQ(y[0], -20.0);
    EXPECT_DOUBLE_EQ(z[0], -30.0);
    
    // Last voxel should be at origin + (dimensions-1) * spacing
    EXPECT_DOUBLE_EQ(x[4], -10.0 + 4 * 2.0);  // -2.0
    EXPECT_DOUBLE_EQ(y[4], -20.0 + 4 * 2.0);  // -12.0
    EXPECT_DOUBLE_EQ(z[4], -30.0 + 4 * 2.0);  // -22.0
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
