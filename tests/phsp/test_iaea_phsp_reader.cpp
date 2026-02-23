#include <gtest/gtest.h>
#include "phsp/IAEAPhspReader.hpp"
#include "phsp/IAEAHeaderParser.hpp"
#include <filesystem>

using namespace optirad;

static std::string getDataDir() {
#ifdef OPTIRAD_DATA_DIR
    return OPTIRAD_DATA_DIR;
#else
    return "data";
#endif
}

class IAEAPhspReaderTest : public ::testing::Test {
protected:
    std::string headerPath;
    std::string phspPath;
    bool filesExist = false;

    void SetUp() override {
        std::string machineDir = getDataDir() + "/machines/Varian_TrueBeam6MV";
        if (std::filesystem::exists(machineDir)) {
            for (auto& entry : std::filesystem::directory_iterator(machineDir)) {
                if (entry.path().extension() == ".IAEAheader") {
                    headerPath = entry.path().string();
                    phspPath = entry.path().stem().string();
                    // Reconstruct full phsp path
                    phspPath = entry.path().parent_path().string() + "/" +
                               entry.path().stem().string() + ".IAEAphsp";
                    if (std::filesystem::exists(phspPath)) {
                        filesExist = true;
                    }
                    break;
                }
            }
        }
    }
};

TEST_F(IAEAPhspReaderTest, ReadSampledParticles) {
    if (!filesExist) {
        GTEST_SKIP() << "No IAEA PSF files found";
    }

    IAEAHeaderParser parser;
    auto info = parser.parse(headerPath);

    IAEAPhspReader reader;
    // Read every 1000th particle
    auto particles = reader.readSampled(phspPath, info, 1000);

    EXPECT_GT(particles.size(), 0) << "Should read some particles";

    // Verify first particle has valid data
    if (!particles.empty()) {
        const auto& p = particles.particles()[0];
        // Position should be in mm (converted from cm)
        // Direction cosines should have magnitude ~1
        double dirMag = p.direction[0] * p.direction[0] +
                        p.direction[1] * p.direction[1] +
                        p.direction[2] * p.direction[2];
        EXPECT_NEAR(dirMag, 1.0, 0.1) << "Direction should be approximately unit vector";
    }
}

TEST_F(IAEAPhspReaderTest, ReadSubset) {
    if (!filesExist) {
        GTEST_SKIP() << "No IAEA PSF files found";
    }

    IAEAHeaderParser parser;
    auto info = parser.parse(headerPath);

    IAEAPhspReader reader;
    auto particles = reader.readSubset(phspPath, info, 0, 100);

    EXPECT_EQ(particles.size(), 100) << "Should read exactly 100 particles";
}

TEST_F(IAEAPhspReaderTest, NonExistentFileThrows) {
    IAEAHeaderParser parser;
    IAEAHeaderInfo info;
    info.recordLength = 25;
    info.totalParticles = 10;

    IAEAPhspReader reader;
    EXPECT_THROW(reader.readAll("/nonexistent.IAEAphsp", info), std::runtime_error);
}
