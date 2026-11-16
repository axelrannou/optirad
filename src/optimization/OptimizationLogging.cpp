#include "OptimizationLogging.hpp"
#include "utils/Logger.hpp"
#include "objectives/DVHObjective.hpp"
#include "objectives/SquaredDeviation.hpp"
#include "objectives/SquaredOverdose.hpp"
#include "geometry/Structure.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace optirad {

bool shouldLogOptimizationProgress(int step, int interval) {
    if (interval <= 0) throw std::invalid_argument("Optimization log interval must be positive");
    return step > 0 && (step == 1 || step % interval == 0);
}

void logDaoProgressHeader(int maxIterations, size_t memorySize) {
    Logger::info("=== DAO Optimizer (Memory: " + std::to_string(memorySize) +
                 ", Max Iter: " + std::to_string(maxIterations) + ") ===");
    std::ostringstream out;
    out << std::setw(5) << "Iter" << ' ' << std::setw(5) << "Acc"
        << ' ' << std::setw(12) << "Objective" << ' ' << std::setw(10) << "GradInf"
        << ' ' << std::setw(10) << "Step" << ' ' << std::setw(4) << "LS"
        << ' ' << std::setw(10) << "MU" << ' ' << std::setw(10) << "Time(s)";
    Logger::info(out.str());
    Logger::info(std::string(73, '-'));
}

void logDaoProgress(const DaoProgressInfo& info) {
    std::ostringstream out;
    out << std::setw(5) << info.iteration << ' ' << std::setw(5) << info.acceptedSteps
        << ' ' << std::scientific << std::setprecision(4) << std::setw(12) << info.objective
        << ' ' << std::setw(10) << info.preStepGradInf
        << ' ' << std::setprecision(2) << std::setw(10) << info.stepSize
        << ' ' << std::setw(4) << info.lsIters
        << ' ' << std::defaultfloat << std::setprecision(5) << std::setw(10) << info.totalMu
        << ' ' << std::setprecision(4) << std::setw(10) << info.elapsedSec;
    Logger::info(out.str());
}

void logObjectiveConfiguration(const std::string& context,
                               const std::vector<ObjectiveFunction*>& objectives) {
    std::ostringstream out;
    out << context << ": " << objectives.size() << " weighted dose objectives (soft penalties)";
    for (size_t i = 0; i < objectives.size(); ++i)
        out << "\n  [" << i + 1 << "] " << objectives[i]->describe();
    Logger::info(out.str());
    for (const auto& warning : objectiveConfigurationWarnings(objectives))
        Logger::warn(context + ": " + warning);
}

std::vector<std::string> objectiveConfigurationWarnings(const std::vector<ObjectiveFunction*>& objectives) {
    std::vector<std::string> warnings;
    std::vector<const Structure*> targets, boundedTargets;
    for (size_t i = 0; i < objectives.size(); ++i) {
        const auto* objective = objectives[i];
        const std::string label = "[" + std::to_string(i + 1) + "] " + objective->describe();
        bool limitsDoseAbove = dynamic_cast<const SquaredOverdose*>(objective) ||
                               dynamic_cast<const SquaredDeviation*>(objective);
        if (const auto* dvh = dynamic_cast<const DVHObjective*>(objective)) {
            const bool isMax = dvh->getType() == DVHObjective::Type::MAX_DVH;
            if (isMax && dvh->getVolumeFraction() >= 1.0)
                warnings.push_back(label + " has no effect: at most 100% of a structure above a dose is always true;"
                                   " use volume 0% for a maximum dose (Eclipse Upper 0%)");
            else if (!isMax && dvh->getVolumeFraction() <= 0.0)
                warnings.push_back(label + " has no effect: at least 0% of a structure above a dose is always true;"
                                   " use volume 100% for a minimum dose (Eclipse Lower 100%)");
            limitsDoseAbove = isMax && dvh->getVolumeFraction() < 1.0;
        }
        const auto* structure = objective->getStructure();
        if (!structure || !structure->isTarget()) continue;
        if (std::find(targets.begin(), targets.end(), structure) == targets.end()) targets.push_back(structure);
        if (limitsDoseAbove) boundedTargets.push_back(structure);
    }
    for (const auto* target : targets)
        if (std::find(boundedTargets.begin(), boundedTargets.end(), target) == boundedTargets.end())
            warnings.push_back("target " + target->getName() + " has no upper objective: nothing limits its"
                               " maximum dose (add MaxDVH volume 0% at ~105% Rx, Eclipse Upper)");
    return warnings;
}

OptimizationBreakdown computeOptimizationBreakdown(
                              const std::vector<ObjectiveFunction*>& objectives,
                              const std::vector<double>& dose,
                              const std::vector<std::pair<std::string, double>>& penalties) {
    OptimizationBreakdown result;
    auto& contributions = result.contributions;
    for (size_t i = 0; i < objectives.size(); ++i)
        contributions.emplace_back("[" + std::to_string(i + 1) + "] " + objectives[i]->describe(),
                                   objectives[i]->compute(dose));
    for (double value : dose) result.maxDose = std::max(result.maxDose, value);
    contributions.insert(contributions.end(), penalties.begin(), penalties.end());
    for (const auto& entry : contributions) result.total += entry.second;
    return result;
}

void logOptimizationBreakdown(const std::string& context, const OptimizationBreakdown& breakdown) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(6) << context << ": weighted cost breakdown total=" << breakdown.total
        << ", max dose=" << breakdown.maxDose << " Gy";
    for (const auto& entry : breakdown.contributions) {
        out << "\n  " << entry.first << ": cost=" << entry.second
            << " share=" << std::fixed << std::setprecision(2)
            << (breakdown.total > 0.0 ? 100.0 * entry.second / breakdown.total : 0.0) << "%"
            << std::scientific << std::setprecision(6);
    }
    Logger::info(out.str());
}

void logOptimizationBreakdown(const std::string& context,
                              const std::vector<ObjectiveFunction*>& objectives,
                              const std::vector<double>& dose,
                              const std::vector<std::pair<std::string, double>>& penalties) {
    logOptimizationBreakdown(context, computeOptimizationBreakdown(objectives, dose, penalties));
}

} // namespace optirad
