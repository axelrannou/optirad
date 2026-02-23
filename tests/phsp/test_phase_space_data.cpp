#include <gtest/gtest.h>
#include "phsp/PhaseSpaceData.hpp"
#include <cmath>

using namespace optirad;

class PhaseSpaceDataTest : public ::testing::Test {
protected:
    PhaseSpaceData createTestData() {
        PhaseSpaceData data;
        // Create 100 test particles
        for (int i = 0; i < 100; ++i) {
            Particle p;
            p.type = (i % 10 == 0) ? ParticleType::Electron : ParticleType::Photon;
            p.position = {static_cast<double>(i) - 50.0, 0.0, 100.0};
            p.direction = {0.0, 0.0, 1.0};
            p.energy = 1.0 + static_cast<double>(i) * 0.05;
            p.weight = 1.0;
            data.addParticle(p);
        }
        return data;
    }
};

TEST_F(PhaseSpaceDataTest, AddAndCount) {
    auto data = createTestData();
    EXPECT_EQ(data.size(), 100);
}

TEST_F(PhaseSpaceDataTest, CountByType) {
    auto data = createTestData();
    auto photons = data.countByType(ParticleType::Photon);
    auto electrons = data.countByType(ParticleType::Electron);
    auto positrons = data.countByType(ParticleType::Positron);
    EXPECT_EQ(photons, 90);
    EXPECT_EQ(electrons, 10);
    EXPECT_EQ(positrons, 0);
}

TEST_F(PhaseSpaceDataTest, ComputeMetrics) {
    auto data = createTestData();
    auto metrics = data.computeMetrics();

    EXPECT_EQ(metrics.totalCount, 100);
    EXPECT_EQ(metrics.photonCount, 90);
    EXPECT_EQ(metrics.electronCount, 10);
    EXPECT_DOUBLE_EQ(metrics.positronCount, 0);

    EXPECT_GT(metrics.meanEnergy, 0.0);
    EXPECT_GE(metrics.maxEnergy, metrics.minEnergy);
}

TEST_F(PhaseSpaceDataTest, Sample) {
    auto data = createTestData();
    auto sampled = data.sample(10, 42);
    EXPECT_EQ(sampled.size(), 10);
}

TEST_F(PhaseSpaceDataTest, SampleLargerThanDataReturnsAll) {
    auto data = createTestData();
    auto sampled = data.sample(200, 42);
    EXPECT_EQ(sampled.size(), 100u);
}

TEST_F(PhaseSpaceDataTest, FilterByType) {
    auto data = createTestData();
    data.filterByType(ParticleType::Electron);
    EXPECT_EQ(data.size(), 10);
    for (const auto& p : data.particles()) {
        EXPECT_EQ(p.type, ParticleType::Electron);
    }
}

TEST_F(PhaseSpaceDataTest, FilterByJaws) {
    // Create particles at Z=27cm scoring plane (like TrueBeam) with positions in mm
    PhaseSpaceData data;
    for (int i = 0; i < 100; ++i) {
        Particle p;
        p.type = ParticleType::Photon;
        p.position = {static_cast<double>(i - 50) * 0.1, 0.0, 271.7}; // mm, X in [-5,4.9]
        p.direction = {0.0, 0.0, 0.999}; // mostly forward
        p.energy = 2.0;
        p.weight = 1.0;
        data.addParticle(p);
    }
    // Jaw convention: jawX1/Y1 = positive opening from center (maps to -jawX1)
    // jawX2/Y2 = positive opening from center (maps to +jawX2)
    // SAD=1000mm, scoring plane at 27.17 cm
    // projScale = 1000 / (27.17 * 10) = 3.68
    // projX = X * 3.68, range [-18.4, 18.0]
    // With jaws: X from -200 to +200, Y from -200 to +200 → all particles pass
    data.filterByJaws(200.0, 200.0, 200.0, 200.0, 1000.0, 27.17);
    EXPECT_GT(data.size(), 0u);
    EXPECT_LE(data.size(), 100u);
}
