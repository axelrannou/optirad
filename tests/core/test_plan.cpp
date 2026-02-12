#include <gtest/gtest.h>
#include "core/Plan.hpp"

namespace optirad::tests {

class PlanTest : public ::testing::Test {
protected:
    Plan plan;
};

TEST_F(PlanTest, InitialNameIsEmpty) {
    EXPECT_EQ(plan.getName(), "");
}

TEST_F(PlanTest, SetNameWorks) {
    plan.setName("Test Plan");
    EXPECT_EQ(plan.getName(), "Test Plan");
}

TEST_F(PlanTest, InitialBeamCountIsZero) {
    EXPECT_EQ(plan.getNumBeams(), 0);
}

TEST_F(PlanTest, AddBeamIncreasesCount) {
    Beam beam;
    beam.setGantryAngle(45.0);
    plan.addBeam(beam);
    EXPECT_EQ(plan.getNumBeams(), 1);
}

TEST_F(PlanTest, GetBeamsWorks) {
    Beam beam1;
    beam1.setGantryAngle(0.0);
    plan.addBeam(beam1);
    
    Beam beam2;
    beam2.setGantryAngle(90.0);
    plan.addBeam(beam2);
    
    auto& beams = plan.getBeams();
    EXPECT_EQ(beams.size(), 2);
    EXPECT_DOUBLE_EQ(beams[0].getGantryAngle(), 0.0);
    EXPECT_DOUBLE_EQ(beams[1].getGantryAngle(), 90.0);
}

TEST_F(PlanTest, SetFractionsWorks) {
    plan.setNumOfFractions(5);
    EXPECT_EQ(plan.getNumOfFractions(), 5);
}

TEST_F(PlanTest, SetRadiationModeWorks) {
    plan.setRadiationMode("6MV");
    EXPECT_EQ(plan.getRadiationMode(), "6MV");
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
