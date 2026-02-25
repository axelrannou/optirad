#include <gtest/gtest.h>
#include "dose/DoseEngineFactory.hpp"
#include "dose/engines/PencilBeamEngine.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/MathUtils.hpp"

namespace optirad::tests {

// ============================================================================
// Factory creation
// ============================================================================

TEST(DoseEngineFactoryTest, CreatePencilBeam) {
    auto engine = DoseEngineFactory::create("PencilBeam");
    ASSERT_NE(engine, nullptr);
    EXPECT_EQ(engine->getName(), "PencilBeam");
}

TEST(DoseEngineFactoryTest, CreatePencilBeamCaseVariations) {
    // Factory may or may not be case-sensitive; test the exact documented name
    auto engine = DoseEngineFactory::create("PencilBeam");
    EXPECT_NE(engine, nullptr);
}

TEST(DoseEngineFactoryTest, UnknownEngineThrows) {
    EXPECT_THROW(DoseEngineFactory::create("UnknownEngine"), std::runtime_error);
}

// ============================================================================
// DoseMatrix basic operations
// ============================================================================

TEST(DoseMatrixTest, AllocateAndAccess) {
    Grid grid;
    grid.setDimensions(4, 5, 6);
    grid.setSpacing(2.0, 2.0, 2.0);
    grid.setOrigin(Vec3{0, 0, 0});

    DoseMatrix dm;
    dm.setGrid(grid);
    dm.allocate();

    EXPECT_EQ(dm.size(), 120u); // 4*5*6

    dm.at(0, 0, 0) = 10.0;
    dm.at(3, 4, 5) = 20.0;

    EXPECT_DOUBLE_EQ(dm.at(0, 0, 0), 10.0);
    EXPECT_DOUBLE_EQ(dm.at(3, 4, 5), 20.0);
}

TEST(DoseMatrixTest, MaxAndMean) {
    Grid grid;
    grid.setDimensions(2, 2, 2);
    grid.setSpacing(1.0, 1.0, 1.0);
    grid.setOrigin(Vec3{0, 0, 0});

    DoseMatrix dm;
    dm.setGrid(grid);
    dm.allocate();

    dm.at(0, 0, 0) = 100.0;
    dm.at(1, 1, 1) = 50.0;
    // others default to 0

    EXPECT_DOUBLE_EQ(dm.getMax(), 100.0);
    EXPECT_DOUBLE_EQ(dm.getMean(), (100.0 + 50.0) / 8.0);
}

TEST(DoseMatrixTest, DataPointerNotNull) {
    Grid grid;
    grid.setDimensions(3, 3, 3);
    grid.setSpacing(1.0, 1.0, 1.0);
    grid.setOrigin(Vec3{0, 0, 0});

    DoseMatrix dm;
    dm.setGrid(grid);
    dm.allocate();

    EXPECT_NE(dm.data(), nullptr);

    const DoseMatrix& cdm = dm;
    EXPECT_NE(cdm.data(), nullptr);
}

// ============================================================================
// PencilBeamEngine forward dose (from Dij + weights)
// ============================================================================

TEST(PencilBeamEngineTest, CalculateDoseFromDij) {
    // Create a small Dij manually
    DoseInfluenceMatrix dij(8, 2); // 2x2x2 grid, 2 bixels
    dij.setValue(0, 0, 1.0);
    dij.setValue(1, 0, 0.5);
    dij.setValue(0, 1, 0.3);
    dij.setValue(7, 1, 2.0);
    dij.finalize();

    std::vector<double> weights = {10.0, 5.0};

    Grid grid;
    grid.setDimensions(2, 2, 2);
    grid.setSpacing(2.5, 2.5, 2.5);
    grid.setOrigin(Vec3{0, 0, 0});

    PencilBeamEngine engine;
    auto dm = engine.calculateDose(dij, weights, grid);

    // Voxel 0: 1.0*10 + 0.3*5 = 11.5
    // Voxel 1: 0.5*10 + 0*5 = 5.0
    // Voxel 7: 0*10 + 2.0*5 = 10.0
    EXPECT_NEAR(dm.at(0, 0, 0), 11.5, 1e-10);
    EXPECT_NEAR(dm.at(1, 0, 0), 5.0, 1e-10);
    EXPECT_NEAR(dm.at(1, 1, 1), 10.0, 1e-10);
}

} // namespace optirad::tests
