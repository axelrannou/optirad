#include <gtest/gtest.h>

#include "optimization/OptimizationLogging.hpp"
#include "optimization/objectives/DVHObjective.hpp"
#include "optimization/objectives/SquaredDeviation.hpp"
#include "optimization/objectives/SquaredOverdose.hpp"
#include "optimization/objectives/SquaredUnderdose.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "dose/DoseInfluenceMatrix.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace optirad::tests {
namespace {

size_t countOccurrences(const std::string& text, const std::string& needle) {
    size_t count = 0, position = 0;
    while ((position = text.find(needle, position)) != std::string::npos) {
        ++count;
        position += needle.size();
    }
    return count;
}

class CountingDoseOperator : public IDoseOperator {
public:
    explicit CountingDoseOperator(const DoseInfluenceMatrix& dij) : m_dij(dij) {}

    size_t getNumVoxels() const override { return m_dij.getNumVoxels(); }
    size_t getNumBixels() const override { return m_dij.getNumBixels(); }
    std::vector<double> computeDose(const std::vector<double>& weights) const override {
        ++doseCalls;
        return m_dij.computeDose(weights);
    }
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                    std::vector<double>& grad) const override {
        m_dij.accumulateTransposeProduct(voxelGrad, grad);
    }
    double getMaxValue() const override { return m_dij.getMaxValue(); }
    bool hasAllRows() const override { return m_dij.hasAllRows(); }
    bool isRowComputed(size_t voxel) const override { return m_dij.isRowComputed(voxel); }

    mutable size_t doseCalls = 0;

private:
    const DoseInfluenceMatrix& m_dij;
};

double loggedValue(const std::string& log, const std::string& field) {
    const auto position = log.find(field);
    if (position == std::string::npos) {
        ADD_FAILURE() << "Missing field: " << field << "\n" << log;
        return 0.0;
    }
    return std::stod(log.substr(position + field.size()));
}

} // namespace

TEST(OptimizationLoggingTest, ProgressCadenceIncludesFirstStepAndEveryConfiguredInterval) {
    for (int step : {1, 50, 100, 150})
        EXPECT_TRUE(shouldLogOptimizationProgress(step, 50)) << step;
    for (int step : {0, 2, 10, 49, 51, 99})
        EXPECT_FALSE(shouldLogOptimizationProgress(step, 50)) << step;
    EXPECT_TRUE(shouldLogOptimizationProgress(10, 10));
    EXPECT_FALSE(shouldLogOptimizationProgress(11, 10));
    EXPECT_THROW(shouldLogOptimizationProgress(1, 0), std::invalid_argument);
}

TEST(OptimizationLoggingTest, DaoProgressTableAlignsColumnsAndKeepsPrecisionAcrossRows) {
    const auto flags = std::cout.flags();
    const auto precision = std::cout.precision();
    DaoProgressInfo info{1, 1, 1.23456789e8, 138.5488, 1e-12, 7, 1743.26, 76.12345};
    testing::internal::CaptureStdout();
    logDaoProgressHeader(500, 10);
    logDaoProgress(info);
    info.iteration = info.acceptedSteps = 50;
    info.objective = 4.987654321e5;
    info.elapsedSec = 120.54321;
    logDaoProgress(info);
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_EQ(std::cout.flags(), flags);
    EXPECT_EQ(std::cout.precision(), precision);

    std::istringstream input(log);
    std::vector<std::string> lines;
    for (std::string line; std::getline(input, line);) lines.push_back(line);
    ASSERT_EQ(lines.size(), 5u);
    EXPECT_EQ(lines[0], "[INFO] === DAO Optimizer (Memory: 10, Max Iter: 500) ===");
    const std::string prefix = "[INFO] ";
    const std::array<size_t, 8> widths{5, 5, 12, 10, 10, 4, 10, 10};
    const std::array<std::string, 8> columns{
        "Iter", "Acc", "Objective", "GradInf", "Step", "LS", "MU", "Time(s)"};
    const std::array<std::string, 8> values{
        "1", "1", "1.2346e+08", "1.3855e+02", "1.00e-12", "7", "1743.3", "76.12"};
    for (size_t i = 1; i < lines.size(); ++i) {
        EXPECT_EQ(lines[i].size(), 80u);
        EXPECT_EQ(lines[i].find(prefix), 0u);
    }
    EXPECT_EQ(lines[2], prefix + std::string(73, '-'));
    size_t position = prefix.size();
    for (size_t i = 0; i < widths.size(); ++i) {
        EXPECT_EQ(lines[1].substr(position, widths[i]),
                  std::string(widths[i] - columns[i].size(), ' ') + columns[i]);
        EXPECT_EQ(lines[3].substr(position, widths[i]),
                  std::string(widths[i] - values[i].size(), ' ') + values[i]);
        position += widths[i] + 1;
    }
    EXPECT_NE(lines[4].find("4.9877e+05"), std::string::npos);
    EXPECT_NE(lines[4].find("1.00e-12"), std::string::npos);
    EXPECT_NE(lines[4].find("120.5"), std::string::npos);
}

TEST(OptimizationLoggingTest, DaoProgressStaysCompactWithLargeValuesAndTinySteps) {
    const int iterations = std::numeric_limits<int>::max();
    const double large = 1.23456789e300;
    testing::internal::CaptureStdout();
    logDaoProgress({iterations, iterations, large, large, std::numeric_limits<double>::denorm_min(),
                    iterations, large, large});
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_LE(log.size() - 1, 100u);
    EXPECT_NE(log.find("4.94e-324"), std::string::npos);
    EXPECT_NE(log.find("1.2346e+300"), std::string::npos);
}

TEST(OptimizationLoggingTest, DescribesEveryObjectiveTypeAndDvhVolume) {
    Structure structure;
    structure.setName("Target");
    structure.setVoxelIndices({0, 1});
    SquaredDeviation deviation;
    deviation.setPrescribedDose(60.0);
    SquaredOverdose overdose;
    overdose.setMaxDose(65.0);
    SquaredUnderdose underdose;
    underdose.setMinDose(57.0);
    DVHObjective minDvh, maxDvh;
    minDvh.setDoseThreshold(59.0);
    minDvh.setVolumeFraction(0.98);
    maxDvh.setType(DVHObjective::Type::MAX_DVH);
    maxDvh.setDoseThreshold(63.0);
    maxDvh.setVolumeFraction(0.02);
    std::vector<ObjectiveFunction*> objectives{&deviation, &overdose, &underdose, &minDvh, &maxDvh};
    for (auto* objective : objectives) {
        objective->setStructure(&structure);
        objective->setWeight(12.0);
        objective->setVoxelIndices({1});
    }
    testing::internal::CaptureStdout();
    logObjectiveConfiguration("Test", objectives);
    const auto log = testing::internal::GetCapturedStdout();
    for (const auto* expected : {
            "5 weighted dose objectives (soft penalties)",
            "Target / SquaredDeviation dose=60 Gy",
            "Target / SquaredOverdose dose=65 Gy",
            "Target / SquaredUnderdose dose=57 Gy",
            "Target / MinDVH dose=59 Gy volume=98%",
            "Target / MaxDVH dose=63 Gy volume=2%",
            "weight=12 voxels=1"})
        EXPECT_NE(log.find(expected), std::string::npos) << expected;
}

TEST(OptimizationLoggingTest, ReportsWeightedContributionsAndExactHotspotCost) {
    SquaredDeviation deviation;
    deviation.setPrescribedDose(1.0);
    deviation.setWeight(2.0);
    deviation.setVoxelIndices({0, 1});
    testing::internal::CaptureStdout();
    logOptimizationBreakdown("Test", {&deviation}, {2.0, 3.0}, 1.0, 1.0, 2.0,
                             {{"Extra", 1.0}});
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(loggedValue(log, "weight=2 voxels=2: cost="), 5.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, "Global hotspot: cost="), 694.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, "Extra: cost="), 1.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, "breakdown total="), 700.0);
    EXPECT_NE(log.find("hotspot voxels=2"), std::string::npos);
    EXPECT_NE(log.find("share=99.14%"), std::string::npos);
}

TEST(OptimizationLoggingTest, ReportsEachMixedObjectiveWithItsOwnWeightAndNormalization) {
    SquaredDeviation deviation;
    deviation.setPrescribedDose(1.0);
    SquaredOverdose overdose;
    overdose.setMaxDose(3.0);
    overdose.setWeight(2.0);
    SquaredUnderdose underdose;
    underdose.setMinDose(3.0);
    underdose.setWeight(4.0);
    DVHObjective minDvh, maxDvh;
    minDvh.setDoseThreshold(3.0);
    minDvh.setVolumeFraction(0.75);
    minDvh.setWeight(3.0);
    maxDvh.setType(DVHObjective::Type::MAX_DVH);
    maxDvh.setDoseThreshold(2.0);
    maxDvh.setVolumeFraction(0.25);
    maxDvh.setWeight(5.0);
    std::vector<ObjectiveFunction*> objectives{&deviation, &overdose, &underdose, &minDvh, &maxDvh};
    for (auto* objective : objectives) objective->setVoxelIndices({0, 1, 2, 3});
    testing::internal::CaptureStdout();
    logOptimizationBreakdown("Test", objectives, {0.0, 2.0, 4.0, 6.0}, 0.0, 0.0, 0.0);
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(loggedValue(log, deviation.describe() + ": cost="), 9.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, overdose.describe() + ": cost="), 5.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, underdose.describe() + ": cost="), 10.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, minDvh.describe() + ": cost="), 3.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, maxDvh.describe() + ": cost="), 20.0);
    EXPECT_DOUBLE_EQ(loggedValue(log, "breakdown total="), 47.0);
}

TEST(OptimizationLoggingTest, ZeroCostAndDisabledHotspotAreReportedWithoutNan) {
    SquaredDeviation deviation;
    deviation.setVoxelIndices({0});
    testing::internal::CaptureStdout();
    logOptimizationBreakdown("Test", {&deviation}, {0.0, 10.0}, 1.0, 1.0, 0.0);
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(loggedValue(log, "breakdown total="), 0.0);
    EXPECT_NE(log.find("share=0.00%"), std::string::npos);
    EXPECT_EQ(log.find("nan"), std::string::npos);
    EXPECT_EQ(log.find("inf"), std::string::npos);
    EXPECT_NE(log.find("hotspot voxels=0"), std::string::npos);
}

TEST(OptimizationLoggingTest, LbfgsBreakdownMatchesTheSolverIncludingRegularizers) {
    DoseInfluenceMatrix dij(2, 2);
    dij.setValue(0, 0, 1.0);
    dij.setValue(1, 1, 1.0);
    dij.finalize();
    SquaredDeviation deviation;
    deviation.setPrescribedDose(1.0);
    deviation.setWeight(2.0);
    deviation.setVoxelIndices({0, 1});
    LBFGSOptimizer optimizer;
    optimizer.setMaxIterations(0);
    optimizer.setInitialWeights({2.0, 3.0});
    optimizer.setPrescriptionDose(1.0);
    optimizer.setHotspotThreshold(1.0);
    optimizer.setHotspotPenalty(2.0);
    optimizer.setSpatialSmoothing(0.5, {{1}, {0}});
    optimizer.setL2Regularization(0.25);
    optimizer.setL1Regularization(0.1);
    testing::internal::CaptureStdout();
    const auto result = optimizer.optimize(dij, {&deviation}, {});
    const auto log = testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(result.finalObjective, 703.25);
    EXPECT_DOUBLE_EQ(loggedValue(log, "breakdown total="), result.finalObjective);
    EXPECT_DOUBLE_EQ(loggedValue(log, "Spatial smoothing: cost="), 0.5);
    EXPECT_DOUBLE_EQ(loggedValue(log, "L2 regularization: cost="), 3.25);
    EXPECT_DOUBLE_EQ(loggedValue(log, "L1 / total weight: cost="), 0.5);
    EXPECT_NE(log.find("L-BFGS-B initial"), std::string::npos);
    EXPECT_NE(log.find("L-BFGS-B final"), std::string::npos);
    EXPECT_EQ(countOccurrences(log, "weighted cost breakdown total="), 2u);
}

TEST(OptimizationLoggingTest, VerboseLoggingPreservesResultsAndReportsAtTheExpectedCadence) {
    DoseInfluenceMatrix dij(2, 2);
    dij.setValue(0, 0, 0.5);
    dij.setValue(1, 1, 1.0);
    dij.finalize();
    CountingDoseOperator countedDij(dij);
    SquaredDeviation deviation;
    deviation.setPrescribedDose(2.0);
    deviation.setVoxelIndices({0, 1});
    LBFGSOptimizer optimizer;
    optimizer.setMaxIterations(20);
    optimizer.setHotspotPenalty(0.0);
    optimizer.setTolerance(0.0);
    std::vector<int> callbacks;
    optimizer.setIterationCallback([&](const IterationInfo& info) { callbacks.push_back(info.iteration); });
    testing::internal::CaptureStdout();
    const auto verbose = optimizer.optimize(countedDij, {&deviation}, {});
    const auto log = testing::internal::GetCapturedStdout();
    const size_t verboseDoseCalls = countedDij.doseCalls;
    ASSERT_GE(callbacks.size(), 2u);
    EXPECT_EQ(callbacks[0], 1);
    EXPECT_EQ(callbacks[1], 10);
    const auto verboseCallbacks = callbacks;
    EXPECT_NE(log.find("\n       1"), std::string::npos);
    EXPECT_NE(log.find("\n      10"), std::string::npos);
    EXPECT_EQ(countOccurrences(log, "weighted cost breakdown total="), 2u);
    EXPECT_EQ(log.find("L-BFGS-B iteration "), std::string::npos);
    optimizer.setVerbose(false);
    countedDij.doseCalls = 0;
    callbacks.clear();
    testing::internal::CaptureStdout();
    const auto quiet = optimizer.optimize(countedDij, {&deviation}, {Constraint{}});
    const auto quietLog = testing::internal::GetCapturedStdout();
    EXPECT_EQ(verbose.weights, quiet.weights);
    EXPECT_DOUBLE_EQ(verbose.finalObjective, quiet.finalObjective);
    EXPECT_EQ(verbose.iterations, quiet.iterations);
    EXPECT_EQ(verboseCallbacks, callbacks);
    EXPECT_EQ(verboseDoseCalls, countedDij.doseCalls + 2);
    EXPECT_EQ(log.find("not enforced"), std::string::npos);
    EXPECT_NE(quietLog.find("Constraint objects supplied but not enforced"), std::string::npos);
    EXPECT_EQ(quietLog.find("weighted cost breakdown"), std::string::npos);
}

} // namespace optirad::tests
