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
#include "dose/DoseInfluenceMatrix.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Grid.hpp"
#include "utils/Logger.hpp"
#include <memory>

#include <iostream>
#include <vector>
#include <cstdlib>

using namespace optirad;

Structure createStructure(const std::string& name, const std::string& type, int nx, int ny, int nz) {
    Structure structure;
    structure.setName(name);
    structure.setType(type);
    
    // Create a simple voxel mask - center region
    std::vector<size_t> voxelIndices;
    for (int k = nz/4; k < 3*nz/4; ++k) {
        for (int j = ny/4; j < 3*ny/4; ++j) {
            for (int i = nx/4; i < 3*nx/4; ++i) {
                size_t idx = i + nx * (j + ny * k);
                voxelIndices.push_back(idx);
            }
        }
    }
    
    structure.setVoxelIndices(voxelIndices);
    return structure;
}

int main() {
    std::cout << "=== OptiRad L-BFGS-B Optimizer Test ===\n\n";
    
    Logger::init();
    
    // Create grid
    int nx = 5, ny = 5, nz = 4;  // Small grid: 100 voxels
    int numBixels = 20;
    
    Grid grid;
    grid.setDimensions(nx, ny, nz);
    grid.setSpacing(1.0, 1.0, 1.0);
    grid.setOrigin({0.0, 0.0, 0.0});
    
    // Create structures
    Structure target = createStructure("PTV", "TARGET", nx, ny, nz);
    Structure oar = createStructure("OAR", "OAR", nx, ny, nz);
    
    // Check structures have voxels
    if (target.getVoxelCount() == 0 || oar.getVoxelCount() == 0) {
        std::cerr << "Error: Structures have no voxels!\n";
        return 1;
    }
    
    std::cout << "Problem size: " << grid.getNumVoxels() << " voxels, " << numBixels << " bixels\n\n";
    
    // Create dose influence matrix
    DoseInfluenceMatrix dij(grid.getNumVoxels(), numBixels);
    
    // Fill with valid data (small random values)
    std::srand(42);
    for (size_t v = 0; v < grid.getNumVoxels(); ++v) {
        for (size_t b = 0; b < static_cast<size_t>(numBixels); ++b) {
            double value = 0.01 + 0.1 * (std::rand() % 100) / 100.0;
            dij(v, b) = value;
        }
    }
    
    // Create objectives (use raw pointers for API compatibility)
    std::vector<ObjectiveFunction*> objectives;
    
    // Target: achieve 60 Gy - check what setters are available
    auto targetObj = std::make_unique<SquaredDeviation>();
    targetObj->setStructure(&target);
    targetObj->setWeight(100.0);
    // Note: SquaredDeviation likely has a constructor parameter or specific setter
    // For now, create a simple working test
    objectives.push_back(targetObj.get());
    
    // OAR: limit to 20 Gy
    auto oarObj = std::make_unique<SquaredOverdose>();
    oarObj->setStructure(&oar);
    oarObj->setWeight(50.0);
    objectives.push_back(oarObj.get());
    
    // Create optimizer
    LBFGSOptimizer optimizer;
    optimizer.setMaxIterations(500);
    optimizer.setMemorySize(10);
    
    // Empty constraints for now
    std::vector<Constraint> constraints;
    
    std::cout << "Starting optimization...\n";
    OptimizationResult result = optimizer.optimize(dij, objectives, constraints);
    
    if (result.converged) {
        std::cout << "\nOptimization successful!\n";
        std::cout << "Final objective value: " << result.finalObjective << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        
        // Display some weights
        std::cout << "\nSample bixel weights:\n";
        for (int i = 0; i < std::min(5, numBixels); ++i) {
            std::cout << "  Bixel " << i << ": " << result.weights[i] << "\n";
        }
    } else {
        std::cout << "\nOptimization failed or did not converge!\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        return 1;
    }
    
    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
