#include "OptimizationLogging.hpp"
#include "utils/Logger.hpp"
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
}

void logOptimizationBreakdown(const std::string& context,
                              const std::vector<ObjectiveFunction*>& objectives,
                              const std::vector<double>& dose,
                              double prescriptionDose, double hotspotThreshold, double hotspotPenalty,
                              const std::vector<std::pair<std::string, double>>& penalties) {
    std::vector<std::pair<std::string, double>> contributions;
    for (size_t i = 0; i < objectives.size(); ++i)
        contributions.emplace_back("[" + std::to_string(i + 1) + "] " + objectives[i]->describe(),
                                   objectives[i]->compute(dose));
    double hotspot = 0.0, maxDose = 0.0;
    size_t hotspotVoxels = 0;
    for (double value : dose) {
        maxDose = std::max(maxDose, value);
        if (prescriptionDose <= 0.0 || hotspotPenalty <= 0.0) continue;
        const double excess = value - prescriptionDose * hotspotThreshold;
        if (excess <= 0.0) continue;
        ++hotspotVoxels;
        const double normalized = excess / prescriptionDose;
        hotspot += hotspotPenalty * normalized * normalized * (1.0 + 10.0 * normalized) * excess * excess;
    }
    contributions.emplace_back("Global hotspot", hotspot);
    contributions.insert(contributions.end(), penalties.begin(), penalties.end());
    double total = 0.0;
    for (const auto& entry : contributions) total += entry.second;

    std::ostringstream out;
    out << std::scientific << std::setprecision(6) << context << ": weighted cost breakdown total=" << total
        << ", max dose=" << maxDose << " Gy, hotspot voxels=" << hotspotVoxels;
    for (const auto& entry : contributions) {
        out << "\n  " << entry.first << ": cost=" << entry.second
            << " share=" << std::fixed << std::setprecision(2)
            << (total > 0.0 ? 100.0 * entry.second / total : 0.0) << "%"
            << std::scientific << std::setprecision(6);
    }
    Logger::info(out.str());
}

} // namespace optirad
