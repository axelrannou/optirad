#pragma once

#include "../IOptimizer.hpp"
#include <deque>
#include <chrono>

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

    // L-BFGS specific parameters
    void setMemorySize(int m);
    void setMaxFluence(double maxFluence);
    void setVerbose(bool verbose);

private:
    // Parameters
    int m_maxIterations = 500;
    double m_tolerance = 1e-5;
    int m_memorySize = 10;
    double m_maxFluence = 10.0;
    bool m_verbose = true;
    int m_progressEvery = 10;
    
    // Line search parameters
    int m_lineSearchMaxIter = 20;
    double m_c1 = 1e-4;      // Armijo condition
    double m_c2 = 0.9;       // Wolfe curvature condition
    double m_initialStepSize = 1.0;
    
    // L-BFGS history
    std::deque<std::vector<double>> m_sHistory;
    std::deque<std::vector<double>> m_yHistory;
    std::deque<double> m_rhoHistory;
    int m_currentMemorySize = 0;
    
    // Compute total objective and gradient
    double computeObjectiveAndGradient(
        const std::vector<double>& weights,
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        std::vector<double>& grad
    );
    
    // L-BFGS two-loop recursion
    void computeLBFGSDirection(const std::vector<double>& grad, std::vector<double>& dir);
    
    // Project direction for bound constraints
    void projectSearchDirection(
        const std::vector<double>& x,
        std::vector<double>& dir,
        double lb, double ub
    );
    
    // Line search with Wolfe conditions
    double lineSearch(
        const std::vector<double>& x,
        double fval,
        const std::vector<double>& grad,
        const std::vector<double>& dir,
        double lb, double ub,
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        std::vector<double>& x_new,
        double& fval_new,
        std::vector<double>& grad_new,
        int& lsIter
    );
    
    // Update L-BFGS history
    void updateHistory(const std::vector<double>& s, const std::vector<double>& y);
    
    // Reset history
    void resetHistory();
    
    // Compute projected gradient norm
    double computeProjectedGradNorm(
        const std::vector<double>& x,
        const std::vector<double>& grad,
        double lb, double ub
    );
    
    // Vector operations
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    static void axpy(double a, const std::vector<double>& x, std::vector<double>& y);
    static void scale(double a, std::vector<double>& x);
};

} // namespace optirad
