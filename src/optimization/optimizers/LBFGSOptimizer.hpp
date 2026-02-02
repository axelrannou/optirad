#pragma once

#include "../IOptimizer.hpp"

namespace optirad {

class LBFGSOptimizer : public IOptimizer {
public:
    std::string getName() const override;

    OptimizationResult optimize(
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        const std::vector<Constraint>& constraints
    ) override;

    void setMaxIterations(int maxIter) override;
    void setTolerance(double tol) override;

private:
    int m_maxIterations = 500;
    double m_tolerance = 1e-6;
};

} // namespace optirad
