#include "LBFGSOptimizer.hpp"
#include "utils/Logger.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

std::string LBFGSOptimizer::getName() const { return "LBFGS"; }

void LBFGSOptimizer::setMaxIterations(int maxIter) { m_maxIterations = maxIter; }
void LBFGSOptimizer::setTolerance(double tol) { m_tolerance = tol; }
void LBFGSOptimizer::setMemorySize(int m) { m_memorySize = m; }
void LBFGSOptimizer::setMaxFluence(double maxFluence) { m_maxFluence = maxFluence; }
void LBFGSOptimizer::setVerbose(bool verbose) { m_verbose = verbose; }

OptimizationResult LBFGSOptimizer::optimize(
    const DoseInfluenceMatrix& dij,
    const std::vector<ObjectiveFunction*>& objectives,
    const std::vector<Constraint>& constraints
) {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    Logger::info("Starting L-BFGS-B optimization...");
    Logger::info("Max iterations: " + std::to_string(m_maxIterations));
    Logger::info("Memory size: " + std::to_string(m_memorySize));
    
    OptimizationResult result;
    result.converged = false;
    result.iterations = 0;
    
    // Initialize
    resetHistory();
    
    int n = static_cast<int>(dij.getNumBixels());
    double lb = 0.0;
    double ub = m_maxFluence;
    
    // Initialize weights
    std::vector<double> x(n, 1.0);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        x[i] = std::max(lb, std::min(ub, x[i]));
    }
    
    // Compute initial objective and gradient
    std::vector<double> grad(n);
    double fval = computeObjectiveAndGradient(x, dij, objectives, grad);
    
    std::vector<double> objectiveHistory;
    objectiveHistory.reserve(m_maxIterations + 1);
    objectiveHistory.push_back(fval);
    
    if (m_verbose) {
        std::cout << "\n=== L-BFGS-B Optimizer (Memory: " << m_memorySize << ") ===\n";
        std::cout << std::setw(8) << "Iter" 
                  << std::setw(16) << "Objective" 
                  << std::setw(16) << "ProjGrad" 
                  << std::setw(12) << "StepSize" 
                  << std::setw(10) << "LS Iters" << "\n";
        std::cout << std::string(62, '-') << "\n";
    }
    
    std::vector<double> dir(n);
    std::vector<double> x_new(n);
    std::vector<double> grad_new(n);
    std::vector<double> s(n);
    std::vector<double> y(n);
    
    // Main optimization loop
    for (int iter = 1; iter <= m_maxIterations; ++iter) {
        result.iterations = iter;
        
        // Compute projected gradient norm
        double projGradNorm = computeProjectedGradNorm(x, grad, lb, ub);
        
        // Check convergence
        if (projGradNorm < m_tolerance) {
            result.converged = true;
            if (m_verbose) {
                std::cout << std::setw(8) << iter 
                          << std::setw(16) << std::scientific << fval 
                          << std::setw(16) << projGradNorm 
                          << "  (Converged)\n";
            }
            break;
        }
        
        // Check relative objective change
        if (iter > 1 && objectiveHistory.size() >= 2) {
            double prev = objectiveHistory[objectiveHistory.size() - 2];
            if (prev != 0) {
                double relObjChange = std::abs(fval - prev) / std::abs(prev);
                if (relObjChange < 1e-7) {
                    result.converged = true;
                    if (m_verbose) {
                        std::cout << std::setw(8) << iter 
                                  << std::setw(16) << std::scientific << fval 
                                  << std::setw(16) << projGradNorm 
                                  << "  (RelObj converged)\n";
                    }
                    break;
                }
            }
        }
        
        // Compute search direction using L-BFGS
        computeLBFGSDirection(grad, dir);
        
        // Project search direction for bound constraints
        projectSearchDirection(x, dir, lb, ub);
        
        // Line search
        int lsIter = 0;
        double alpha = lineSearch(x, fval, grad, dir, lb, ub, dij, objectives, 
                                  x_new, fval, grad_new, lsIter);
        
        // Display progress
        if (m_verbose && (iter % m_progressEvery == 0 || iter == 1)) {
            double initialObj = objectiveHistory[0];
            double improvement = 0;
            if (initialObj > 0) {
                improvement = 100.0 * (initialObj - fval) / initialObj;
            }
            std::cout << std::setw(8) << iter 
                      << std::setw(16) << std::scientific << fval 
                      << std::setw(16) << projGradNorm 
                      << std::setw(12) << std::fixed << alpha 
                      << std::setw(10) << lsIter 
                      << "  (" << std::fixed << std::setprecision(1) << improvement << "% imp)\n";
        }
        
        // Compute s and y for L-BFGS update
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            s[i] = x_new[i] - x[i];
            y[i] = grad_new[i] - grad[i];
        }
        
        // Check curvature condition
        double sDotY = dot(s, y);
        double sNorm = std::sqrt(dot(s, s));
        double yNorm = std::sqrt(dot(y, y));
        
        if (sDotY > 1e-10 * sNorm * yNorm) {
            updateHistory(s, y);
        }
        
        // Update current point
        std::swap(x, x_new);
        std::swap(grad, grad_new);
        
        objectiveHistory.push_back(fval);
        
        // Adaptive restart if progress stalls
        if (iter > 10) {
            size_t histSize = objectiveHistory.size();
            double recentObj = objectiveHistory[std::max(size_t(0), histSize - 6)];
            double recentProgress = std::abs(recentObj - fval) / (std::abs(recentObj) + 1e-10);
            if (recentProgress < 1e-6) {
                resetHistory();
                if (m_verbose) {
                    std::cout << "         (L-BFGS history reset due to slow progress)\n";
                }
            }
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(endTime - startTime).count();
    
    result.weights = std::move(x);
    result.finalObjective = fval;
    
    if (m_verbose) {
        std::cout << "\nOptimization completed in " << std::fixed << std::setprecision(2) 
                  << elapsed << "s after " << result.iterations << " iterations\n";
        std::cout << "Final objective: " << std::scientific << result.finalObjective << "\n";
        std::cout << "Status: " << (result.converged ? "Converged" : "Max iterations") << "\n";
    }
    
    Logger::info("Optimization complete.");
    return result;
}

double LBFGSOptimizer::computeObjectiveAndGradient(
    const std::vector<double>& weights,
    const DoseInfluenceMatrix& dij,
    const std::vector<ObjectiveFunction*>& objectives,
    std::vector<double>& grad
) {
    int numVoxels = static_cast<int>(dij.getNumVoxels());
    int n = static_cast<int>(weights.size());
    
    // Compute dose = Dij * weights
    std::vector<double> dose = dij.computeDose(weights);
    
    // Initialize gradient to zero
    std::fill(grad.begin(), grad.end(), 0.0);
    
    // Compute voxel-space gradient
    std::vector<double> gVox(numVoxels, 0.0);
    
    double objVal = 0.0;
    
    // Sum objectives
    for (const auto* obj : objectives) {
        objVal += obj->compute(dose);
        
        std::vector<double> objGrad = obj->gradient(dose);
        #pragma omp parallel for
        for (int i = 0; i < numVoxels; ++i) {
            gVox[i] += objGrad[i];
        }
    }
    
    // Compute gradient using optimized transpose product: grad = Dij^T * gVox
    dij.accumulateTransposeProduct(gVox, grad);
    
    return objVal;
}

void LBFGSOptimizer::computeLBFGSDirection(const std::vector<double>& grad, std::vector<double>& dir) {
    int n = static_cast<int>(grad.size());
    
    // Copy gradient to q
    std::vector<double> q = grad;
    std::vector<double> alphas(m_currentMemorySize);
    
    // First loop (backward)
    for (int i = m_currentMemorySize - 1; i >= 0; --i) {
        alphas[i] = m_rhoHistory[i] * dot(m_sHistory[i], q);
        axpy(-alphas[i], m_yHistory[i], q);
    }
    
    // Initial Hessian approximation (scaled identity)
    if (m_currentMemorySize > 0) {
        double gamma = dot(m_sHistory.back(), m_yHistory.back()) / 
                      dot(m_yHistory.back(), m_yHistory.back());
        scale(gamma, q);
    }
    
    // Second loop (forward)
    for (int i = 0; i < m_currentMemorySize; ++i) {
        double beta = m_rhoHistory[i] * dot(m_yHistory[i], q);
        axpy(alphas[i] - beta, m_sHistory[i], q);
    }
    
    // Search direction is negative gradient direction
    dir.resize(n);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        dir[i] = -q[i];
    }
}

void LBFGSOptimizer::projectSearchDirection(
    const std::vector<double>& x,
    std::vector<double>& dir,
    double lb, double ub
) {
    int n = static_cast<int>(x.size());
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        // At lower bound and direction points outside -> zero it
        if (x[i] <= lb + 1e-10 && dir[i] < 0) {
            dir[i] = 0;
        }
        // At upper bound and direction points outside -> zero it
        else if (x[i] >= ub - 1e-10 && dir[i] > 0) {
            dir[i] = 0;
        }
    }
}

double LBFGSOptimizer::lineSearch(
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
) {
    int n = static_cast<int>(x.size());
    double alpha = m_initialStepSize;
    double alphaLow = 0;
    double alphaHigh = std::numeric_limits<double>::infinity();
    
    double dirDeriv = dot(grad, dir);
    
    // Not a descent direction
    if (dirDeriv >= 0) {
        x_new = x;
        fval_new = fval;
        grad_new = grad;
        lsIter = 0;
        return 0;
    }
    
    for (lsIter = 1; lsIter <= m_lineSearchMaxIter; ++lsIter) {
        // Compute trial point with projection
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            x_new[i] = std::max(lb, std::min(ub, x[i] + alpha * dir[i]));
        }
        
        // Evaluate objective and gradient at trial point
        fval_new = computeObjectiveAndGradient(x_new, dij, objectives, grad_new);
        
        // Check Armijo condition (sufficient decrease)
        if (fval_new > fval + m_c1 * alpha * dirDeriv) {
            alphaHigh = alpha;
            alpha = 0.5 * (alphaLow + alphaHigh);
            continue;
        }
        
        // Check Wolfe curvature condition
        double dirDerivNew = dot(grad_new, dir);
        
        // Strong Wolfe condition
        if (std::abs(dirDerivNew) <= -m_c2 * dirDeriv) {
            break;  // Found acceptable step
        }
        
        if (dirDerivNew < 0) {
            alphaLow = alpha;
            if (std::isinf(alphaHigh)) {
                alpha = 2 * alpha;
            } else {
                alpha = 0.5 * (alphaLow + alphaHigh);
            }
        } else {
            alphaHigh = alpha;
            alpha = 0.5 * (alphaLow + alphaHigh);
        }
    }
    
    return alpha;
}

void LBFGSOptimizer::updateHistory(const std::vector<double>& s, const std::vector<double>& y) {
    double sDotY = dot(y, s);
    
    // Check for zero curvature to prevent NaN propagation
    if (std::abs(sDotY) < 1e-14) {
        Logger::warn("Skipping BFGS update - zero curvature");
        return;
    }
    
    double rho = 1.0 / sDotY;
    
    if (m_currentMemorySize < m_memorySize) {
        m_sHistory.push_back(s);
        m_yHistory.push_back(y);
        m_rhoHistory.push_back(rho);
        m_currentMemorySize++;
    } else {
        // Remove oldest, add newest
        m_sHistory.pop_front();
        m_yHistory.pop_front();
        m_rhoHistory.pop_front();
        m_sHistory.push_back(s);
        m_yHistory.push_back(y);
        m_rhoHistory.push_back(rho);
    }
}

void LBFGSOptimizer::resetHistory() {
    m_sHistory.clear();
    m_yHistory.clear();
    m_rhoHistory.clear();
    m_currentMemorySize = 0;
}

double LBFGSOptimizer::computeProjectedGradNorm(
    const std::vector<double>& x,
    const std::vector<double>& grad,
    double lb, double ub
) {
    int n = static_cast<int>(x.size());
    double maxVal = 0.0;
    
    #pragma omp parallel for reduction(max:maxVal)
    for (int i = 0; i < n; ++i) {
        double proj = x[i] - grad[i];
        proj = std::max(lb, std::min(ub, proj));
        double diff = std::abs(x[i] - proj);
        if (diff > maxVal) maxVal = diff;
    }
    
    return maxVal;
}

double LBFGSOptimizer::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    int n = static_cast<int>(a.size());
    
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < n; ++i) {
        result += a[i] * b[i];
    }
    
    return result;
}

void LBFGSOptimizer::axpy(double a, const std::vector<double>& x, std::vector<double>& y) {
    int n = static_cast<int>(x.size());
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        y[i] += a * x[i];
    }
}

void LBFGSOptimizer::scale(double a, std::vector<double>& x) {
    int n = static_cast<int>(x.size());
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        x[i] *= a;
    }
}

} // namespace optirad
