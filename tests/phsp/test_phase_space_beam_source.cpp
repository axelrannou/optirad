#include <gtest/gtest.h>
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "core/Machine.hpp"
#include <filesystem>

using namespace optirad;

static std::string getDataDir() {
#ifdef OPTIRAD_DATA_DIR
    return OPTIRAD_DATA_DIR;
#else
    return "data";
#endif
}

class PhaseSpaceBeamSourceTest : public ::testing::Test {
protected:
    bool machineAvailable = false;
    Machine machine;

    void SetUp() override {
        std::string machineDir = getDataDir() + "/machines/Varian_TrueBeam6MV";
        if (!std::filesystem::exists(machineDir)) return;

        // Check if PSF files exist
        bool hasHeader = false, hasPhsp = false;
        for (auto& entry : std::filesystem::directory_iterator(machineDir)) {
            if (entry.path().extension() == ".IAEAheader") hasHeader = true;
            if (entry.path().extension() == ".IAEAphsp") hasPhsp = true;
        }
        if (!hasHeader || !hasPhsp) return;

        // Construct a minimal phase-space machine
        MachineMeta meta;
        meta.name = "Varian_TrueBeam6MV";
        meta.SAD = 1000.0;
        meta.machineType = MachineType::PhaseSpace;
        machine.setMeta(meta);

        MachineGeometry geom;
        geom.jawX1Min = -200.0;
        geom.jawX1Max = 0.0;
        geom.jawX2Min = 0.0;
        geom.jawX2Max = 200.0;
        geom.jawY1Min = -200.0;
        geom.jawY1Max = 0.0;
        geom.jawY2Min = 0.0;
        geom.jawY2Max = 200.0;
        geom.beamEnergyMV = 6.0;
        geom.phaseSpaceDir = machineDir;

        // Discover PSF files
        for (auto& entry : std::filesystem::directory_iterator(machineDir)) {
            if (entry.path().extension() == ".IAEAheader") {
                geom.phaseSpaceFileNames.push_back(entry.path().stem().string());
            }
        }
        geom.numPhaseSpaceFiles = static_cast<int>(geom.phaseSpaceFileNames.size());

        machine.setGeometry(geom);
        machineAvailable = true;
    }
};

TEST_F(PhaseSpaceBeamSourceTest, ConfigureAndBuild) {
    if (!machineAvailable) {
        GTEST_SKIP() << "Varian TrueBeam 6MV PSF files not available";
    }

    PhaseSpaceBeamSource source;
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    source.configure(machine, 0.0, 0.0, 0.0, iso);

    // Build with a small number of particles for speed
    source.build(10000, 1000);

    EXPECT_TRUE(source.isBuilt());

    const auto& metrics = source.getMetrics();
    EXPECT_GT(metrics.totalCount, 0);
    EXPECT_GT(metrics.photonCount, 0);

    const auto& vizSample = source.getVisualizationSample();
    EXPECT_GT(vizSample.size(), 0);
    EXPECT_LE(vizSample.size(), 1000u);
}

TEST_F(PhaseSpaceBeamSourceTest, EnergyHistogram) {
    if (!machineAvailable) {
        GTEST_SKIP() << "Varian TrueBeam 6MV PSF files not available";
    }

    PhaseSpaceBeamSource source;
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    source.configure(machine, 0.0, 0.0, 0.0, iso);
    source.build(10000, 1000);

    auto hist = source.computeEnergyHistogram(20);
    EXPECT_EQ(hist.size(), 20);

    // At least one bin should have counts
    int64_t totalCounts = 0;
    for (const auto& [energy, count] : hist) {
        totalCounts += count;
    }
    EXPECT_GT(totalCounts, 0);
}

TEST_F(PhaseSpaceBeamSourceTest, BuildNotCalledIsNotBuilt) {
    PhaseSpaceBeamSource source;
    EXPECT_FALSE(source.isBuilt());
}
