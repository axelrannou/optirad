#include <gtest/gtest.h>

#include "optimization/objectives/NormalTissueObjective.hpp"
#include "geometry/Structure.hpp"

#include <cmath>
#include <memory>

namespace optirad::tests {

namespace {

Grid makeGrid(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz) {
    Grid grid;
    grid.setDimensions(nx, ny, nz);
    grid.setSpacing(dx, dy, dz);
    grid.setOrigin({0.0, 0.0, 0.0});
    return grid;
}

void addStructure(StructureSet& set, const std::string& name, const std::string& type, std::vector<size_t> voxels) {
    auto structure = std::make_unique<Structure>();
    structure->setName(name);
    structure->setType(type);
    structure->setVoxelIndices(voxels);
    set.addStructure(std::move(structure));
}

NtoSettings settings() {
    NtoSettings nto;
    nto.priority = 3.0;
    nto.distanceFromTargetMm = 4.0;
    nto.startDosePct = 100.0;
    nto.endDosePct = 50.0;
    nto.fallOffPerMm = 0.1;
    return nto;
}

} // namespace

TEST(NormalTissueObjectiveTest, LimitStaysAtTheStartDoseThenFallsOffTowardsTheEndDose) {
    const auto nto = settings();
    EXPECT_DOUBLE_EQ(NormalTissueObjective::limit(nto, 50.0, 0.0), 50.0);
    EXPECT_DOUBLE_EQ(NormalTissueObjective::limit(nto, 50.0, 4.0), 50.0);
    const double decay = std::exp(-1.0);
    EXPECT_NEAR(NormalTissueObjective::limit(nto, 50.0, 14.0), 50.0 * decay + 25.0 * (1.0 - decay), 1e-12);
    EXPECT_NEAR(NormalTissueObjective::limit(nto, 50.0, 1e4), 25.0, 1e-9);
}

TEST(NormalTissueObjectiveTest, CoversTheBodyOutsideTheTargetsWithDistanceBasedLimits) {
    const Grid grid = makeGrid(11, 1, 1, 2.0, 2.0, 2.0);
    StructureSet structures;
    addStructure(structures, "PTV", "PTV", {0, 1});
    addStructure(structures, "BODY", "EXTERNAL", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    const auto nto = NormalTissueObjective::build(settings(), 50.0, &structures, grid, grid);
    ASSERT_NE(nto, nullptr);
    const auto& voxels = nto->voxels();
    ASSERT_EQ(voxels.size(), 8u);   // body minus target; voxel 10 is outside the body
    for (size_t k = 0; k < voxels.size(); ++k) {
        EXPECT_EQ(voxels[k], k + 2);
        const double distance = 2.0 * static_cast<double>(voxels[k] - 1) - 1.0;   // from the border, not the centre
        EXPECT_NEAR(nto->limits()[k], NormalTissueObjective::limit(settings(), 50.0, distance), 1e-9);
    }
    EXPECT_DOUBLE_EQ(nto->getWeight(), 3.0);
    EXPECT_NE(nto->describe().find("NTO ("), std::string::npos);
}

TEST(NormalTissueObjectiveTest, DistanceUsesTheVoxelSpacingOfEveryAxis) {
    const Grid grid = makeGrid(5, 5, 1, 1.0, 2.0, 3.0);
    StructureSet structures;
    addStructure(structures, "PTV", "PTV", {2 + 5 * 2});
    NtoSettings nto = settings();
    nto.distanceFromTargetMm = 0.0;
    const auto objective = NormalTissueObjective::build(nto, 50.0, &structures, grid, grid);
    ASSERT_NE(objective, nullptr);   // no body: the whole grid outside the target
    ASSERT_EQ(objective->voxels().size(), 24u);
    auto limitAt = [&](size_t x, size_t y) {
        const size_t voxel = x + 5 * y;
        const auto& voxels = objective->voxels();
        for (size_t k = 0; k < voxels.size(); ++k)
            if (voxels[k] == voxel) return objective->limits()[k];
        return -1.0;
    };
    EXPECT_NEAR(limitAt(4, 2), NormalTissueObjective::limit(nto, 50.0, 1.5), 1e-9);
    EXPECT_NEAR(limitAt(2, 4), NormalTissueObjective::limit(nto, 50.0, 3.5), 1e-9);
    EXPECT_NEAR(limitAt(4, 4), NormalTissueObjective::limit(nto, 50.0, std::sqrt(20.0) - 0.5), 1e-9);
}

TEST(NormalTissueObjectiveTest, DisabledOrWithoutTargetBuildsNothing) {
    const Grid grid = makeGrid(4, 1, 1, 1.0, 1.0, 1.0);
    StructureSet structures;
    addStructure(structures, "BODY", "EXTERNAL", {0, 1, 2, 3});
    EXPECT_EQ(NormalTissueObjective::build(settings(), 50.0, &structures, grid, grid), nullptr);
    addStructure(structures, "PTV", "PTV", {0});
    NtoSettings off = settings();
    off.enabled = false;
    EXPECT_EQ(NormalTissueObjective::build(off, 50.0, &structures, grid, grid), nullptr);
    EXPECT_EQ(NormalTissueObjective::build(settings(), 0.0, &structures, grid, grid), nullptr);
    EXPECT_NE(NormalTissueObjective::build(settings(), 50.0, &structures, grid, grid), nullptr);
}

TEST(NormalTissueObjectiveTest, CostPenalizesOnlyDoseAboveTheLimitAndGradientMatchesFiniteDifferences) {
    NormalTissueObjective nto({0, 2, 3}, {10.0, 20.0, 5.0});
    nto.setWeight(2.0);
    std::vector<double> dose = {12.0, 100.0, 15.0, 9.0};
    EXPECT_NEAR(nto.compute(dose), 2.0 * (4.0 + 0.0 + 16.0), 1e-12);   // sum, like MaxDVH 0%
    const auto gradient = nto.gradient(dose);
    ASSERT_EQ(gradient.size(), dose.size());
    EXPECT_DOUBLE_EQ(gradient[1], 0.0);   // not a normal tissue voxel
    EXPECT_DOUBLE_EQ(gradient[2], 0.0);   // below its limit
    const double h = 1e-6;
    for (size_t voxel = 0; voxel < dose.size(); ++voxel) {
        auto up = dose, down = dose;
        up[voxel] += h;
        down[voxel] -= h;
        EXPECT_NEAR(gradient[voxel], (nto.compute(up) - nto.compute(down)) / (2.0 * h), 1e-6) << voxel;
    }
    DoseObjectiveSpec description;
    ASSERT_TRUE(nto.gpuDescription(description));
    EXPECT_EQ(description.type, 1);
    EXPECT_DOUBLE_EQ(description.weight / description.indices.size(), 2.0);   // the GPU averages type 1
    EXPECT_EQ(description.thresholds, nto.limits());
}

} // namespace optirad::tests
