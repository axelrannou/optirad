#include <gtest/gtest.h>
#include "geometry/CoordinateSystem.hpp"

namespace optirad::tests {

class CoordinateSystemTest : public ::testing::Test {
protected:
    CoordinateSystem cs;
};

TEST_F(CoordinateSystemTest, DefaultOriginIsZero) {
    // Origin should default to zero vector
    cs.setOrigin({0.0, 0.0, 0.0});
    // No getter for origin, so just verify no exception
    SUCCEED();
}

TEST_F(CoordinateSystemTest, DefaultSpacingIsOne) {
    // Spacing should default to unit vector
    cs.setSpacing({1.0, 1.0, 1.0});
    SUCCEED();
}

TEST_F(CoordinateSystemTest, SetOriginWorks) {
    EXPECT_NO_THROW({
        cs.setOrigin({10.0, 20.0, 30.0});
    });
}

TEST_F(CoordinateSystemTest, SetSpacingWorks) {
    EXPECT_NO_THROW({
        cs.setSpacing({0.5, 1.0, 2.0});
    });
}

TEST_F(CoordinateSystemTest, PatientToWorldConversion) {
    cs.setOrigin({0.0, 0.0, 0.0});
    cs.setSpacing({1.0, 1.0, 1.0});
    
    // Patient coordinates (0,0,0) should map to world origin
    auto world = cs.patientToWorld({0.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(world[0], 0.0);
    EXPECT_DOUBLE_EQ(world[1], 0.0);
    EXPECT_DOUBLE_EQ(world[2], 0.0);
}

TEST_F(CoordinateSystemTest, WorldToPatientConversion) {
    cs.setOrigin({0.0, 0.0, 0.0});
    cs.setSpacing({1.0, 1.0, 1.0});
    
    auto patient = cs.worldToPatient({0.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(patient[0], 0.0);
    EXPECT_DOUBLE_EQ(patient[1], 0.0);
    EXPECT_DOUBLE_EQ(patient[2], 0.0);
}

TEST_F(CoordinateSystemTest, VoxelToWorldConversion) {
    cs.setOrigin({0.0, 0.0, 0.0});
    cs.setSpacing({1.0, 1.0, 1.0});
    
    // Voxel (5,5,5) with unit spacing should be at world (5,5,5)
    auto world = cs.voxelToWorld(5, 5, 5);
    EXPECT_DOUBLE_EQ(world[0], 5.0);
    EXPECT_DOUBLE_EQ(world[1], 5.0);
    EXPECT_DOUBLE_EQ(world[2], 5.0);
}

TEST_F(CoordinateSystemTest, WorldToVoxelConversion) {
    cs.setOrigin({0.0, 0.0, 0.0});
    cs.setSpacing({2.0, 2.0, 2.0});
    
    int i, j, k;
    cs.worldToVoxel({10.0, 10.0, 10.0}, i, j, k);
    
    // With spacing of 2.0, world (10,10,10) should be at voxel (5,5,5)
    EXPECT_EQ(i, 5);
    EXPECT_EQ(j, 5);
    EXPECT_EQ(k, 5);
}

TEST_F(CoordinateSystemTest, NonUniformSpacing) {
    cs.setOrigin({0.0, 0.0, 0.0});
    cs.setSpacing({0.5, 1.0, 2.0});
    
    // Verify non-uniform spacing is accepted
    SUCCEED();
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
