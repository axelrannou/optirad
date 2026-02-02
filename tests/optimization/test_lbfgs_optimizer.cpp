/**
 * test_lbfgs_optimizer.cpp
 * 
 * Integration test for the L-BFGS-B optimizer.
 * Creates a synthetic TPS optimization problem and verifies the optimizer works.
 * 
 * Compile standalone:
 *   g++ -O3 -fopenmp -std=c++17 -I../../src test_lbfgs_optimizer.cpp \
 *       ../../src/optimization/optimizers/LBFGSOptimizer.cpp \
 *       ../../src/optimization/ObjectiveFunction.cpp \
 *       ../../src/optimization/objectives/SquaredDeviation.cpp \
 *       ../../src/optimization/objectives/SquaredOverdose.cpp \
 *       ../../src/optimization/objectives/SquaredUnderdose.cpp \
 *       ../../src/dose/DoseInfluenceMatrix.cpp \
 *       ../../src/geometry/Structure.cpp \
 *       ../../src/utils/Logger.cpp \
 *       -o test_lbfgs_optimizer
 */

#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/objectives/SquaredDeviation.hpp"
#include "optimization/objectives/SquaredOverdose.hpp"
#include "optimization/objectives/SquaredUnderdose.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "geometry/Structure.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace optirad;

/**
 * Create a synthetic dose influence matrix (Dij)
 * Each voxel receives dose from ~3 nearby beamlets
 */
DoseInfluenceMatrix createSyntheticDij(int numVoxels, int numBixels) {
    DoseInfluenceMatrix dij;
    dij.setDimensions(numVoxels, numBixels);
    
    // Each voxel receives dose from a few beamlets
    for (int i = 0; i < numVoxels; ++i) {
        int bixel1 = i % numBixels;
        int bixel2 = (i + 1) % numBixels;
        int bixel3 = (i + 2) % numBixels;
        
        dij.setValue(i, bixel1, 1.0);
        dij.setValue(i, bixel2, 0.5);
        dij.setValue(i, bixel3, 0.3);
    }
    
    return dij;
}

/**
 * Create structure with voxel mask
 */
Structure createStructure(const std::string& name, const std::string& type,
                          int startVoxel, int endVoxel, int totalVoxels) {
    Structure structure;
    structure.setName(name);
    structure.setType(type);
    
    std::vector<bool> mask(totalVoxels, false);
    for (int i = startVoxel; i < endVoxel; ++i) {
        mask[i] = true;
    }
    structure.setMask(mask);
    
    return structure;
}

int main() {
    std::cout << "=== OptiRad L-BFGS-B Optimizer Test ===\n\n";
    
    Logger::init();
    
    // Problem dimensions
    const int numVoxels = 100;
    const int numBixels = 20;
    
    std::cout << "Problem size: " << numVoxels << " voxels, " << numBixels << " bixels\n\n";
    
    // Create synthetic Dij matrix
    DoseInfluenceMatrix dij = createSyntheticDij(numVoxels, numBixels);
    
    // Create structures
    // Target: voxels 0-29 (should receive 2.0 Gy)
    Structure target = createStructure("PTV", "TARGET", 0, 30, numVoxels);
    
    // OAR: voxels 50-69 (should receive < 1.0 Gy)
    Structure oar = createStructure("Rectum", "OAR", 50, 70, numVoxels);
    
    // Create objectives
    SquaredDeviation targetObj;
    targetObj.setStructure(&target);
    targetObj.setPrescribedDose(2.0);
    targetObj.setWeight(100.0);  // High priority for target
    
    SquaredOverdose oarObj;
    oarObj.setStructure(&oar);
    oarObj.setMaxDose(1.0);
    oarObj.setWeight(50.0);
    
    std::vector<ObjectiveFunction*> objectives = {&targetObj, &oarObj};
    std::vector<Constraint> constraints;  // Empty for now
    
    // Create and configure optimizer
    LBFGSOptimizer optimizer;
    optimizer.setMaxIterations(500);
    optimizer.setTolerance(1e-5);
    optimizer.setMemorySize(10);
    optimizer.setMaxFluence(10.0);
    optimizer.setVerbose(true);
    
    // Run optimization
    std::cout << "Starting optimization...\n";
    OptimizationResult result = optimizer.optimize(dij, objectives, constraints);
    
    // Analyze results
    std::cout << "\n=== Results ===\n";
    std::cout << "Status: " << (result.converged ? "Converged" : "Max iterations") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final objective: " << std::scientific << result.finalObjective << "\n";
    
    // Compute final dose distribution
    std::vector<double> finalDose = dij.computeDose(result.weights);
    
    // Analyze target dose
    double targetMean = 0, targetMin = 1e9, targetMax = 0;
    auto targetIndices = target.getVoxelIndices();
    for (size_t idx : targetIndices) {
        targetMean += finalDose[idx];
        targetMin = std::min(targetMin, finalDose[idx]);
        targetMax = std::max(targetMax, finalDose[idx]);
    }
    targetMean /= targetIndices.size();
    
    std::cout << "\nTarget (PTV, prescribed=2.0 Gy):\n";
    std::cout << "  Min:  " << std::fixed << std::setprecision(3) << targetMin << " Gy\n";
    std::cout << "  Mean: " << targetMean << " Gy\n";
    std::cout << "  Max:  " << targetMax << " Gy\n";
    
    // Analyze OAR dose
    double oarMean = 0, oarMax = 0;
    auto oarIndices = oar.getVoxelIndices();
    for (size_t idx : oarIndices) {
        oarMean += finalDose[idx];
        oarMax = std::max(oarMax, finalDose[idx]);
    }
    oarMean /= oarIndices.size();
    
    std::cout << "\nOAR (Rectum, max=1.0 Gy):\n";
    std::cout << "  Mean: " << oarMean << " Gy\n";
    std::cout << "  Max:  " << oarMax << " Gy\n";
    
    // Print optimized weights
    std::cout << "\nOptimized beamlet weights:\n";
    for (int i = 0; i < numBixels; ++i) {
        std::cout << "  w[" << std::setw(2) << i << "] = " 
                  << std::fixed << std::setprecision(4) << result.weights[i] << "\n";
    }
    
    // Basic validation
    bool success = true;
    
    if (!result.converged && result.iterations >= 500) {
        std::cout << "\n[WARNING] Optimizer did not converge within max iterations\n";
    }
    
    if (targetMean < 1.5 || targetMean > 2.5) {
        std::cout << "[FAIL] Target mean dose is too far from prescription\n";
        success = false;
    }
    
    if (oarMax > 2.0) {  // Allow some slack in this simple test
        std::cout << "[FAIL] OAR max dose is too high\n";
        success = false;
    }
    
    std::cout << "\n=== Test " << (success ? "PASSED" : "FAILED") << " ===\n";
    
    return success ? 0 : 1;
}
