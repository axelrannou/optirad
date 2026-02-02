#pragma once

#include "DoseInfluenceMatrix.hpp"
#include "ObjectiveFunction.hpp"
#include "Constraint.hpp"
#include <vector>
#include <string>

namespace optirad {

struct OptimizationResult {
    std::vector<double> weights;
    double finalObjective;
    int iterations;
    bool converged;
};

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
};

} // namespace optirad
