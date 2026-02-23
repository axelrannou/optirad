#include <gtest/gtest.h>
#include "io/MachineLoader.hpp"
#include "core/Machine.hpp"
#include <stdexcept>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad::tests {

class MachineLoaderTest : public ::testing::Test {};

// ============================================================================
// Successful loading
// ============================================================================

TEST_F(MachineLoaderTest, LoadGenericPhotonFromDataDir) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    EXPECT_EQ(machine.getName(), "Generic");
    EXPECT_EQ(machine.getRadiationMode(), "photons");
}

TEST_F(MachineLoaderTest, MetaFieldsCorrect) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& meta = machine.getMeta();
    EXPECT_EQ(meta.name, "Generic");
    EXPECT_EQ(meta.radiationMode, "photons");
    EXPECT_EQ(meta.dataType, "-");
    EXPECT_EQ(meta.createdOn, "27-Oct-2015");
    EXPECT_EQ(meta.createdBy, "wieserh");
    EXPECT_DOUBLE_EQ(meta.SAD, 1000.0);
    EXPECT_DOUBLE_EQ(meta.SCD, 500.0);
}

TEST_F(MachineLoaderTest, DataScalarFieldsCorrect) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& data = machine.getData();
    EXPECT_DOUBLE_EQ(data.energy, 6.0);
    EXPECT_DOUBLE_EQ(data.m, 0.005066);
    EXPECT_DOUBLE_EQ(data.penumbraFWHMatIso, 5.0);
}

TEST_F(MachineLoaderTest, BetasCorrect) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& betas = machine.getData().betas;
    EXPECT_DOUBLE_EQ(betas[0], 0.3252);
    EXPECT_DOUBLE_EQ(betas[1], 0.016);
    EXPECT_DOUBLE_EQ(betas[2], 0.0051);
}

TEST_F(MachineLoaderTest, PrimaryFluenceLoaded) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& fluence = machine.getData().primaryFluence;
    EXPECT_EQ(fluence.size(), 38u);

    // First entry: [0, 1]
    EXPECT_DOUBLE_EQ(fluence[0][0], 0.0);
    EXPECT_DOUBLE_EQ(fluence[0][1], 1.0);

    // Last entry: [325.27, 0.018089]
    EXPECT_NEAR(fluence.back()[0], 325.27, 1e-6);
    EXPECT_NEAR(fluence.back()[1], 0.018089, 1e-6);
}

TEST_F(MachineLoaderTest, KernelPositionsLoaded) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& kp = machine.getData().kernelPos;
    EXPECT_EQ(kp.size(), 360u);
    EXPECT_DOUBLE_EQ(kp[0], 0.0);
    EXPECT_DOUBLE_EQ(kp[1], 0.5);
    EXPECT_DOUBLE_EQ(kp.back(), 179.5);
}

TEST_F(MachineLoaderTest, KernelEntriesLoaded) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& kernel = machine.getData().kernel;

    // Should have 501 entries (SSD 500..1000)
    EXPECT_EQ(kernel.size(), 501u);

    // First entry: SSD = 500
    EXPECT_DOUBLE_EQ(kernel[0].SSD, 500.0);
    EXPECT_EQ(kernel[0].kernel1.size(), 360u);
    EXPECT_EQ(kernel[0].kernel2.size(), 360u);
    EXPECT_EQ(kernel[0].kernel3.size(), 360u);

    // Last entry: SSD = 1000
    EXPECT_DOUBLE_EQ(kernel.back().SSD, 1000.0);
}

TEST_F(MachineLoaderTest, ConstraintsCorrect) {
    Machine machine = MachineLoader::load("photons", "Generic", OPTIRAD_DATA_DIR);
    const auto& c = machine.getConstraints();
    EXPECT_DOUBLE_EQ(c.gantryRotationSpeed[0], 0.0);
    EXPECT_DOUBLE_EQ(c.gantryRotationSpeed[1], 6.0);
    EXPECT_DOUBLE_EQ(c.leafSpeed[0], 0.0);
    EXPECT_DOUBLE_EQ(c.leafSpeed[1], 60.0);
    EXPECT_DOUBLE_EQ(c.monitorUnitRate[0], 1.25);
    EXPECT_DOUBLE_EQ(c.monitorUnitRate[1], 10.0);
}

// ============================================================================
// Error cases
// ============================================================================

TEST_F(MachineLoaderTest, ThrowsOnNonExistentFile) {
    EXPECT_THROW(
        MachineLoader::load("photons", "NonExistent", OPTIRAD_DATA_DIR),
        std::runtime_error
    );
}

TEST_F(MachineLoaderTest, ThrowsOnBadPath) {
    EXPECT_THROW(
        MachineLoader::loadFromFile("/nonexistent/path/machine.json"),
        std::runtime_error
    );
}

// ============================================================================
// Convenience overload (uses compiled-in data dir)
// ============================================================================

TEST_F(MachineLoaderTest, DefaultDataDirOverloadWorks) {
    Machine machine = MachineLoader::load("photons", "Generic");
    EXPECT_EQ(machine.getName(), "Generic");
    EXPECT_DOUBLE_EQ(machine.getSAD(), 1000.0);
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
