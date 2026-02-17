#include <gtest/gtest.h>
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"

namespace optirad::tests {

/**
 * Tests for column-major indexing
 */
class IndexingTest : public ::testing::Test {
protected:
    Grid grid;
    optirad::CTVolume volume;
    
    void SetUp() override {
        // Create a small test grid: 3x4x2 (ny=3, nx=4, nz=2)
        grid.setDimensions(3, 4, 2);  // [ny, nx, nz]
        grid.setSpacing(1.0, 1.0, 1.0);
        grid.setOrigin({0.0, 0.0, 0.0});
        
        volume.setGrid(grid);
        volume.allocate();
        
        // Fill with sequential values for testing
        for (size_t idx = 0; idx < volume.size(); ++idx) {
            volume.data()[idx] = static_cast<int16_t>(idx);
        }
    }
};

TEST_F(IndexingTest, DimensionsAreCorrect) {
    auto dims = grid.getDimensions();
    EXPECT_EQ(dims[0], 3);  // ny
    EXPECT_EQ(dims[1], 4);  // nx
    EXPECT_EQ(dims[2], 2);  // nz
}

TEST_F(IndexingTest, ColumnMajorIndexing) {
    // Test column-major indexing: index = i + j*ny + k*ny*nx
    // With ny=3, nx=4, nz=2
    
    // Voxel (0,0,0) should be at index 0
    EXPECT_EQ(volume.at(0, 0, 0), 0);
    
    // Voxel (1,0,0) should be at index 1 (next row, same column)
    EXPECT_EQ(volume.at(1, 0, 0), 1);
    
    // Voxel (2,0,0) should be at index 2
    EXPECT_EQ(volume.at(2, 0, 0), 2);
    
    // Voxel (0,1,0) should be at index 3 (first row, next column)
    EXPECT_EQ(volume.at(0, 1, 0), 3);
    
    // Voxel (0,0,1) should be at index 12 (first voxel of second slice)
    EXPECT_EQ(volume.at(0, 0, 1), 12);
}

TEST_F(IndexingTest, LinearIndexFormula) {
    size_t ny = 3, nx = 4;
    
    // Test all voxels in first slice
    for (size_t j = 0; j < nx; ++j) {
        for (size_t i = 0; i < ny; ++i) {
            size_t expected_index = i + j * ny;
            int16_t value = volume.at(i, j, 0);
            EXPECT_EQ(value, expected_index) 
                << "at(" << i << "," << j << ",0) should be " << expected_index;
        }
    }
    
    // Test all voxels in second slice
    for (size_t j = 0; j < nx; ++j) {
        for (size_t i = 0; i < ny; ++i) {
            size_t expected_index = i + j * ny + 1 * ny * nx;
            int16_t value = volume.at(i, j, 1);
            EXPECT_EQ(value, expected_index)
                << "at(" << i << "," << j << ",1) should be " << expected_index;
        }
    }
}

TEST_F(IndexingTest, CoordinateArraySizes) {
    auto x = grid.getXCoordinates();
    auto y = grid.getYCoordinates();
    auto z = grid.getZCoordinates();
    
    EXPECT_EQ(x.size(), 4);  // nx
    EXPECT_EQ(y.size(), 3);  // ny
    EXPECT_EQ(z.size(), 2);  // nz
}

TEST_F(IndexingTest, TotalVoxelCount) {
    // Total voxels = ny * nx * nz = 3 * 4 * 2 = 24
    EXPECT_EQ(volume.size(), 24);
    EXPECT_EQ(grid.getNumVoxels(), 24);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
