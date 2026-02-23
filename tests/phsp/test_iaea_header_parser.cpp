#include <gtest/gtest.h>
#include "phsp/IAEAHeaderParser.hpp"
#include <fstream>
#include <filesystem>

using namespace optirad;

static std::string getDataDir() {
#ifdef OPTIRAD_DATA_DIR
    return OPTIRAD_DATA_DIR;
#else
    return "data";
#endif
}

class IAEAHeaderParserTest : public ::testing::Test {
protected:
    std::string headerPath;
    bool headerExists = false;

    void SetUp() override {
        // Look for Varian TrueBeam 6MV IAEA header
        std::string machineDir = getDataDir() + "/machines/Varian_TrueBeam6MV";
        if (std::filesystem::exists(machineDir)) {
            for (auto& entry : std::filesystem::directory_iterator(machineDir)) {
                if (entry.path().extension() == ".IAEAheader") {
                    headerPath = entry.path().string();
                    headerExists = true;
                    break;
                }
            }
        }
    }
};

TEST_F(IAEAHeaderParserTest, ParseHeaderIfAvailable) {
    if (!headerExists) {
        GTEST_SKIP() << "No IAEA header file found in data/machines/Varian_TrueBeam6MV/";
    }

    IAEAHeaderParser parser;
    auto info = parser.parse(headerPath);

    // Basic sanity checks
    EXPECT_GT(info.totalParticles, 0) << "Should have particles";
    EXPECT_GT(info.numPhotons, 0) << "Should have photons";
    EXPECT_GT(info.recordLength, 0) << "Record length should be positive";
    EXPECT_FALSE(info.needsByteSwap()) << "Same endianness expected on x86";
}

TEST_F(IAEAHeaderParserTest, ComputeRecordLength) {
    if (!headerExists) {
        GTEST_SKIP() << "No IAEA header file found";
    }

    IAEAHeaderParser parser;
    auto info = parser.parse(headerPath);

    int computed = info.computeRecordLength();
    EXPECT_EQ(computed, info.recordLength) << "Computed record length should match declared";
}

TEST_F(IAEAHeaderParserTest, NonExistentFileThrows) {
    IAEAHeaderParser parser;
    EXPECT_THROW(parser.parse("/nonexistent/path.IAEAheader"), std::runtime_error);
}

TEST_F(IAEAHeaderParserTest, EnergyStatistics) {
    if (!headerExists) {
        GTEST_SKIP() << "No IAEA header file found";
    }

    IAEAHeaderParser parser;
    auto info = parser.parse(headerPath);

    EXPECT_GT(info.photonMeanEnergy, 0.0) << "Mean photon energy should be positive";
    EXPECT_LE(info.photonMeanEnergy, 10.0) << "Mean photon energy should be < 10 MeV for 6MV";
}
