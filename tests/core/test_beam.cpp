#include <gtest/gtest.h>
#include "core/Beam.hpp"
#include <cmath>

namespace optirad::tests {

class BeamTest : public ::testing::Test {
protected:
    Beam beam;
};

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
    beam.setGantryAngle(450.0);  // Equivalent to 90 degrees
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 450.0);
}

TEST_F(BeamTest, SetIsocenterWorks) {
    beam.setIsocenter(100.0, 200.0, 300.0);
    // Note: There's no getter for isocenter in the header, so we just verify no exception
    SUCCEED();
}

TEST_F(BeamTest, MultipleAngleUpdates) {
    beam.setGantryAngle(0.0);
    beam.setGantryAngle(180.0);
    beam.setGantryAngle(360.0);
    EXPECT_DOUBLE_EQ(beam.getGantryAngle(), 360.0);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
