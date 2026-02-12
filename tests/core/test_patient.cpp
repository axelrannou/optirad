#include <gtest/gtest.h>
#include "core/Patient.hpp"

namespace optirad::tests {

class PatientTest : public ::testing::Test {
protected:
    Patient patient;
};

TEST_F(PatientTest, InitialNameIsEmpty) {
    EXPECT_EQ(patient.getName(), "");
}

TEST_F(PatientTest, InitialIDIsEmpty) {
    EXPECT_EQ(patient.getID(), "");
}

TEST_F(PatientTest, SetNameWorks) {
    patient.setName("John Doe");
    EXPECT_EQ(patient.getName(), "John Doe");
}

TEST_F(PatientTest, SetIDWorks) {
    patient.setID("12345");
    EXPECT_EQ(patient.getID(), "12345");
}

TEST_F(PatientTest, NameCanBeUpdated) {
    patient.setName("Initial");
    EXPECT_EQ(patient.getName(), "Initial");
    
    patient.setName("Updated");
    EXPECT_EQ(patient.getName(), "Updated");
}

TEST_F(PatientTest, IDCanBeUpdated) {
    patient.setID("ID001");
    patient.setID("ID002");
    EXPECT_EQ(patient.getID(), "ID002");
}

TEST_F(PatientTest, EmptyNameCanBeSet) {
    patient.setName("SomeName");
    patient.setName("");
    EXPECT_EQ(patient.getName(), "");
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
