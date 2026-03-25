#pragma once

#include "DoseInfluenceMatrix.hpp"
#include "ObjectiveFunction.hpp"
#include "Constraint.hpp"
#include <vector>
#include <string>
#include <functional>

namespace optirad {

struct OptimizationResult {
    std::vector<double> weights;
    double finalObjective;
    int iterations;
    bool converged;
};

/// Per-iteration progress information reported by the optimizer.
struct IterationInfo {
    int iteration = 0;
    double objective = 0.0;
    double projGradNorm = 0.0;
    double stepSize = 0.0;
    int lsIters = 0;
    double improvement = 0.0;  // % improvement from initial objective
};

using IterationCallback = std::function<void(const IterationInfo&)>;

class IOptimizer {
public:
    virtual ~IOptimizer() = default;

    virtual std::string getName() const = 0;

    virtual OptimizationResult optimize(
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        const std::vector<Constraint>& constraints
    ) = 0;

    virtual void setMaxIterations(int maxIter) = 0;
    virtual void setTolerance(double tol) = 0;

    /// Set an optional callback invoked at regular iteration intervals.
    virtual void setIterationCallback(IterationCallback /*cb*/) {}
};

} // namespace optirad
