#include <gtest/gtest.h>
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"

namespace optirad::tests {

class CTVolumeTest : public ::testing::Test {
protected:
    Grid grid;
    optirad::CTVolume volume;
    
    void SetUp() override {
        grid.setDimensions(10, 10, 10);
        grid.setSpacing(1.0, 1.0, 1.0);
        grid.setOrigin({0.0, 0.0, 0.0});
        
        volume.setGrid(grid);
        volume.allocate();
    }
};

TEST_F(CTVolumeTest, GridIsSet) {
    auto& retrievedGrid = volume.getGrid();
    auto dims = retrievedGrid.getDimensions();
    EXPECT_EQ(dims[0], 10);
    EXPECT_EQ(dims[1], 10);
    EXPECT_EQ(dims[2], 10);
}

TEST_F(CTVolumeTest, AllocateCreatesStorage) {
    size_t expectedSize = 10 * 10 * 10;
    EXPECT_EQ(volume.size(), expectedSize);
}

TEST_F(CTVolumeTest, DataPointerIsValid) {
    auto* ptr = volume.data();
    EXPECT_NE(ptr, nullptr);
}

TEST_F(CTVolumeTest, CanWriteAndReadData) {
    int16_t value = 100;
    volume.at(0, 0, 0) = value;
    EXPECT_EQ(volume.at(0, 0, 0), value);
}

TEST_F(CTVolumeTest, BoundsCanBeAccessedSafely) {
    // Try accessing corners
    EXPECT_NO_THROW({
        volume.at(0, 0, 0);
        volume.at(9, 9, 9);
        volume.at(5, 5, 5);
    });
}

TEST_F(CTVolumeTest, DefaultValuesAreZero) {
    // All values should be zero-initialized
    for (size_t k = 0; k < 10; ++k) {
        for (size_t j = 0; j < 10; ++j) {
            for (size_t i = 0; i < 10; ++i) {
                EXPECT_EQ(volume.at(i, j, k), 0);
            }
        }
    }
}

class DoseVolumeTest : public ::testing::Test {
protected:
    Grid grid;
    optirad::DoseVolume volume;
    
    void SetUp() override {
        grid.setDimensions(5, 5, 5);
        grid.setSpacing(2.0, 2.0, 2.0);
        volume.setGrid(grid);
        volume.allocate();
    }
};

TEST_F(DoseVolumeTest, DoseVolumeWorks) {
    double doseValue = 5.5;
    volume.at(0, 0, 0) = doseValue;
    EXPECT_DOUBLE_EQ(volume.at(0, 0, 0), doseValue);
}

TEST_F(DoseVolumeTest, DoseVolumeSizeIsCorrect) {
    size_t expectedSize = 5 * 5 * 5;
    EXPECT_EQ(volume.size(), expectedSize);
}

TEST_F(DoseVolumeTest, VolumeCoordinateArrays) {
    // Volume should provide coordinate arrays through its grid
    auto x = volume.getXCoordinates();
    auto y = volume.getYCoordinates();
    auto z = volume.getZCoordinates();
    
    EXPECT_EQ(x.size(), 5);
    EXPECT_EQ(y.size(), 5);
    EXPECT_EQ(z.size(), 5);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
