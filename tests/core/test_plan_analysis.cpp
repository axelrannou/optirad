#include "core/workflow/PlanAnalysis.hpp"
#include "gui/AppState.hpp"
#include "gui/panels/DoseStatsCsv.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>

using namespace optirad;

namespace {

Grid makeGrid(size_t nz, double originZ = 0.0, double spacingZ = 1.0) {
    Grid grid;
    grid.setDimensions(1, 1, nz);
    grid.setSpacing(1.0, 1.0, spacingZ);
    grid.setOrigin({0.0, 0.0, originZ});
    return grid;
}

PatientData makePatient(const Grid& grid) {
    PatientData patient;
    auto ct = std::make_unique<Volume<int16_t>>();
    ct->setGrid(grid);
    ct->allocate();
    patient.setCTVolume(std::move(ct));
    auto structures = std::make_unique<StructureSet>();
    auto organ = std::make_unique<Structure>();
    organ->setName("PTV");
    organ->setType("TARGET");
    std::vector<size_t> indices(grid.getNumVoxels());
    std::iota(indices.begin(), indices.end(), 0);
    organ->setVoxelIndices(indices);
    structures->addStructure(std::move(organ));
    patient.setStructureSet(std::move(structures));
    return patient;
}

std::shared_ptr<DoseMatrix> makeDose(const Grid& grid, const std::vector<double>& values) {
    auto dose = std::make_shared<DoseMatrix>();
    dose->setGrid(grid);
    dose->allocate();
    std::copy(values.begin(), values.end(), dose->data());
    return dose;
}

} // namespace

TEST(PlanAnalysis, CommonCoverageUsesIdenticalAnatomyAndRecomputesAllMetrics) {
    auto full = makeGrid(5);
    auto cropped = makeGrid(3, 1.0);
    auto patient = makePatient(full);
    auto fullDose = makeDose(full, {0, 20, 40, 60, 100});
    auto croppedDose = makeDose(cropped, {20, 40, 60});
    const std::vector<Grid> grids{full, cropped};
    auto fullStats = PlanAnalysis::computeStats(*fullDose, patient, full, 40, grids);
    auto croppedStats = PlanAnalysis::computeStats(*croppedDose, patient, cropped, 40, grids);
    ASSERT_EQ(fullStats.size(), 1u);
    ASSERT_EQ(croppedStats.size(), 1u);
    const auto& s = fullStats.front();
    EXPECT_EQ(s.sortedDoses, croppedStats.front().sortedDoses);
    EXPECT_EQ(s.sortedDoses, (std::vector<double>{20, 40, 60}));
    EXPECT_EQ(s.numVoxels, 3u);
    EXPECT_EQ(s.totalStructureVoxels, 5u);
    EXPECT_DOUBLE_EQ(s.coveragePercent(), 60);
    EXPECT_DOUBLE_EQ(s.minDose, 20);
    EXPECT_DOUBLE_EQ(s.maxDose, 60);
    EXPECT_DOUBLE_EQ(s.meanDose, 40);
    EXPECT_NEAR(s.stdDose, std::sqrt(800.0 / 3), 1e-10);
    EXPECT_NEAR(s.d2, 59.2, 1e-10);
    EXPECT_NEAR(s.d5, 58, 1e-10);
    EXPECT_DOUBLE_EQ(s.d50, 40);
    EXPECT_NEAR(s.d95, 22, 1e-10);
    EXPECT_NEAR(s.d98, 20.8, 1e-10);
    EXPECT_DOUBLE_EQ(s.v20, 100);
    EXPECT_NEAR(s.v40, 200.0 / 3, 1e-10);
    EXPECT_NEAR(s.v50, 100.0 / 3, 1e-10);
    EXPECT_NEAR(s.v60, 100.0 / 3, 1e-10);
    EXPECT_NEAR(s.conformityIndex, 2.0 / 3, 1e-10);
    EXPECT_NEAR(s.homogeneityIndex, 0.96, 1e-10);
    const auto curves = PlanAnalysis::computeDVHCurves(fullStats, 60, 4);
    ASSERT_EQ(curves.size(), 1u);
    EXPECT_NEAR(curves.front().volumes[2], 200.0 / 3, 1e-5);
}

TEST(PlanAnalysis, DifferentResolutionInterpolatesOnCTCentresAndRetainsZeros) {
    auto full = makeGrid(5);
    auto coarse = makeGrid(2, 1, 2);
    auto patient = makePatient(full);
    auto dose = makeDose(coarse, {0, 40});
    auto stats = PlanAnalysis::computeStats(*dose, patient, coarse, 60, {full, coarse});
    ASSERT_EQ(stats.size(), 1u);
    // Half-voxel edge support is preserved, rather than extrapolating missing dose.
    EXPECT_EQ(stats.front().sortedDoses, (std::vector<double>{0, 0, 20, 40, 40}));
    EXPECT_EQ(stats.front().numVoxels, 5u);
    EXPECT_DOUBLE_EQ(stats.front().coveragePercent(), 100);
}

TEST(PlanAnalysis, ShiftedAndReversedGridUsesPhysicalIntersection) {
    auto full = makeGrid(5);
    auto reversed = makeGrid(3, 3);
    reversed.setImageOrientation({1, 0, 0, 0, -1, 0});
    auto patient = makePatient(full);
    auto fullDose = makeDose(full, {0, 20, 40, 60, 100});
    auto reversedDose = makeDose(reversed, {60, 40, 20});
    auto a = PlanAnalysis::computeStats(*fullDose, patient, full, 60, {full, reversed});
    auto b = PlanAnalysis::computeStats(*reversedDose, patient, reversed, 60, {full, reversed});
    ASSERT_EQ(a.size(), 1u);
    ASSERT_EQ(b.size(), 1u);
    EXPECT_EQ(a.front().sortedDoses, (std::vector<double>{20, 40, 60}));
    EXPECT_EQ(a.front().sortedDoses, b.front().sortedDoses);
}

TEST(PlanAnalysis, NoCommonCoverageOmitsStatisticsRatherThanReportingZeroDose) {
    auto full = makeGrid(5);
    auto disjoint = makeGrid(2, 10);
    auto patient = makePatient(full);
    auto dose = makeDose(full, {0, 20, 40, 60, 100});
    auto stats = PlanAnalysis::computeStats(*dose, patient, full, 60, {full, disjoint});
    EXPECT_TRUE(stats.empty());
    EXPECT_TRUE(PlanAnalysis::computeDVHCurves(stats, 100).empty());
}

TEST(PlanAnalysis, SingleDoseKeepsExistingSampling) {
    auto full = makeGrid(5);
    auto patient = makePatient(full);
    auto dose = makeDose(full, {0, 20, 40, 60, 100});
    auto stats = PlanAnalysis::computeStats(*dose, patient, full);
    ASSERT_EQ(stats.size(), 1u);
    EXPECT_EQ(stats.front().numVoxels, 5u);
    EXPECT_DOUBLE_EQ(stats.front().meanDose, 44);
    EXPECT_DOUBLE_EQ(stats.front().coveragePercent(), 100);
}

TEST(PlanAnalysis, GuiCacheRefreshesWhenLoadedDoseSetChanges) {
    auto full = makeGrid(5);
    auto cropped = makeGrid(3, 1);
    auto patient = makePatient(full);
    GuiAppState state;
    state.doseStore.addDose("Full", makeDose(full, {0, 20, 40, 60, 100}),
                            std::make_shared<Grid>(full));
    ASSERT_EQ(state.getOrComputeStats(0, patient).size(), 1u);
    EXPECT_EQ(state.getOrComputeStats(0, patient).front().numVoxels, 5u);
    state.doseStore.addDose("Imported", makeDose(cropped, {20, 40, 60}),
                            std::make_shared<Grid>(cropped));
    const auto a = state.getOrComputeStats(0, patient);
    const auto b = state.getOrComputeStats(1, patient);
    EXPECT_EQ(a.front().numVoxels, 3u);
    EXPECT_EQ(a.front().sortedDoses, b.front().sortedDoses);
    auto narrow = makeGrid(1, 2);
    state.doseStore.addDose("Third", makeDose(narrow, {40}), std::make_shared<Grid>(narrow));
    EXPECT_EQ(state.getOrComputeStats(0, patient).front().numVoxels, 1u);
    EXPECT_EQ(state.getOrComputeStats(1, patient).front().numVoxels, 1u);
    EXPECT_EQ(state.getOrComputeStats(2, patient).front().numVoxels, 1u);
    state.doseStore.removeDose(2);
    EXPECT_EQ(state.getOrComputeStats(0, patient).front().numVoxels, 3u);
    state.doseStore.removeDose(1);
    EXPECT_EQ(state.getOrComputeStats(0, patient).front().numVoxels, 5u);
    EXPECT_DOUBLE_EQ(state.getOrComputeStats(0, patient, 20).front().conformityIndex, 0.8);
    EXPECT_DOUBLE_EQ(state.getOrComputeStats(0, patient, 60).front().conformityIndex, 0.4);
}

TEST(PlanAnalysis, CsvExportsSampleCountAndFullStructureCoverage) {
    auto full = makeGrid(5);
    auto cropped = makeGrid(3, 1);
    auto patient = makePatient(full);
    auto dose = makeDose(full, {0, 20, 40, 60, 100});
    auto stats = PlanAnalysis::computeStats(*dose, patient, full, 60, {cropped});
    auto path = std::filesystem::path(testing::TempDir()) / "optirad_plan_analysis_stats.csv";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out);
        writeDoseStatsCsvHeader(out);
        writeDoseStatsCsvRows(out, "Full", stats);
    }
    std::ifstream in(path);
    std::string header, row;
    std::getline(in, header);
    std::getline(in, row);
    EXPECT_NE(header.find("TotalStructureVoxels,Coverage(%)"), std::string::npos);
    EXPECT_EQ(row.substr(0, 18), "Full,PTV,TARGET,3,");
    EXPECT_EQ(row.substr(row.size() - 5), ",5,60");
    in.close();
    std::filesystem::remove(path);
}
