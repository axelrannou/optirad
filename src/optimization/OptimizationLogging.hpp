#pragma once

#include "ObjectiveFunction.hpp"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace optirad {

bool shouldLogOptimizationProgress(int step, int interval);

struct DaoProgressInfo {
    int iteration = 0;
    int acceptedSteps = 0;
    double objective = 0.0;
    double preStepGradInf = 0.0;
    double stepSize = 0.0;
    int lsIters = 0;
    double totalMu = 0.0;
    double elapsedSec = 0.0;
};

void logDaoProgressHeader(int maxIterations, size_t memorySize);
void logDaoProgress(const DaoProgressInfo& info);

void logObjectiveConfiguration(const std::string& context,
                               const std::vector<ObjectiveFunction*>& objectives);

/// Objectives that cannot have any effect, and targets whose maximum dose is not limited.
std::vector<std::string> objectiveConfigurationWarnings(const std::vector<ObjectiveFunction*>& objectives);

struct OptimizationBreakdown {
    std::vector<std::pair<std::string, double>> contributions;
    double total = 0.0;
    double maxDose = 0.0;
};

OptimizationBreakdown computeOptimizationBreakdown(
    const std::vector<ObjectiveFunction*>& objectives, const std::vector<double>& dose,
    const std::vector<std::pair<std::string, double>>& penalties = {});

void logOptimizationBreakdown(const std::string& context, const OptimizationBreakdown& breakdown);

// Diagnostic snapshots only: no gradients or changes to the optimizer state.
void logOptimizationBreakdown(const std::string& context,
                              const std::vector<ObjectiveFunction*>& objectives,
                              const std::vector<double>& dose,
                                                        const std::vector<std::pair<std::string, double>>& penalties = {});

} // namespace optirad
