#include "Plan.hpp"
#include "PatientData.hpp"
#include "geometry/Volume.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"

#include <iostream>
#include <numeric>
#include <vector>
#include <set>

namespace optirad {

void Plan::setName(const std::string& name) { m_name = name; }
const std::string& Plan::getName() const { return m_name; }

void Plan::addBeam(const Beam& beam) { m_beams.push_back(beam); }
const std::vector<Beam>& Plan::getBeams() const { return m_beams; }
size_t Plan::getNumBeams() const { return m_beams.size(); }

std::array<double, 3> Plan::computeIsoCenter() const {
    if (!m_patientData) {
        return {0.0, 0.0, 0.0};
    }

    auto* ct = m_patientData->getCTVolume();
    auto* structureSet = m_patientData->getStructureSet();
    
    if (!ct || !structureSet) {
        return {0.0, 0.0, 0.0};
    }

    const auto& grid = ct->getGrid();
    auto dims = grid.getDimensions();
    auto spacing = grid.getSpacing();
    auto origin = grid.getOrigin();

    // Collect all voxel indices from target structures
    std::set<size_t> targetVoxels;
    int targetCount = 0;
    
    for (size_t i = 0; i < structureSet->getCount(); ++i) {
        const auto* structure = structureSet->getStructure(i);
        if (structure && structure->getType() == "TARGET") {
            targetCount++;
            // Get voxel indices for this target structure
            const auto& voxels = structure->getVoxelIndices();
            if (!voxels.empty()) {
                targetVoxels.insert(voxels.begin(), voxels.end());
            } else {
                std::cerr << "Warning: Target structure '" << structure->getName() 
                          << "' has no voxel indices (contours not rasterized yet)\n";
            }
        }
    }

    if (targetCount == 0) {
        std::cerr << "Warning: No TARGET structures found in structure set\n";
        std::cerr << "Available structures:\n";
        for (size_t i = 0; i < structureSet->getCount(); ++i) {
            const auto* s = structureSet->getStructure(i);
            if (s) {
                std::cerr << "  - " << s->getName() << " (type: " << s->getType() << ")\n";
            }
        }
    }

    if (targetVoxels.empty()) {
        // Fallback: use CT volume center if no target voxels found
        std::cerr << "Warning: No target voxel indices available, using CT volume center as isocenter\n";
        std::cerr << "Note: You may need to rasterize structure contours first\n";
        std::array<double, 3> iso;
        iso[0] = origin[0] + (dims[0] - 1) * spacing[0] / 2.0;
        iso[1] = origin[1] + (dims[1] - 1) * spacing[1] / 2.0;
        iso[2] = origin[2] + (dims[2] - 1) * spacing[2] / 2.0;
        return iso;
    }

    // Convert voxel indices to world coordinates and compute center of gravity
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    
    for (size_t voxelIdx : targetVoxels) {
        // Convert linear index to 3D coordinates (i, j, k)
        size_t k = voxelIdx / (dims[0] * dims[1]);
        size_t j = (voxelIdx % (dims[0] * dims[1])) / dims[0];
        size_t i = voxelIdx % dims[0];
        
        // Convert to world coordinates [mm]
        double x = origin[0] + i * spacing[0];
        double y = origin[1] + j * spacing[1];
        double z = origin[2] + k * spacing[2];
        
        sumX += x;
        sumY += y;
        sumZ += z;
    }

    size_t numVoxels = targetVoxels.size();
    std::array<double, 3> isoCenter;
    isoCenter[0] = sumX / numVoxels;
    isoCenter[1] = sumY / numVoxels;
    isoCenter[2] = sumZ / numVoxels;

    std::cout << "Computed isocenter from " << numVoxels << " target voxels across " 
              << targetCount << " target structure(s)\n";

    return isoCenter;
}

void Plan::printSummary() const {
    std::cout << "\n=== Treatment Plan Summary ===\n";
    std::cout << "  Plan Name:       " << (m_name.empty() ? "(unnamed)" : m_name) << "\n";
    std::cout << "  Radiation Mode:  " << m_radiationMode << "\n";
    std::cout << "  Machine:         " << m_machine.getName() << "\n";
    std::cout << "  SAD:             " << m_machine.getSAD() << " mm\n";
    std::cout << "  SCD:             " << m_machine.getSCD() << " mm\n";
    std::cout << "  Energy:          " << m_machine.getData().energy << " MV\n";
    std::cout << "  Num Fractions:   " << m_numOfFractions << "\n";

    std::cout << "\n  STF Properties:\n";
    std::cout << "    Num Beams:     " << m_stfProperties.numOfBeams << "\n";
    std::cout << "    Bixel Width:   " << m_stfProperties.bixelWidth << " mm\n";

    if (!m_stfProperties.gantryAngles.empty()) {
        std::cout << "    Gantry Angles: [" << m_stfProperties.gantryAngles.front()
                  << " : ... : " << m_stfProperties.gantryAngles.back() << "] deg"
                  << " (" << m_stfProperties.gantryAngles.size() << " angles)\n";

        if (m_stfProperties.gantryAngles.size() > 1) {
            double step = m_stfProperties.gantryAngles[1] - m_stfProperties.gantryAngles[0];
            std::cout << "    Gantry Step:   " << step << " deg\n";
        }
    }

    if (!m_stfProperties.isoCenters.empty()) {
        const auto& iso = m_stfProperties.isoCenters[0];
        std::cout << "    Isocenter:     [" << iso[0] << ", " << iso[1] << ", " << iso[2] << "] mm\n";
    }

    // Machine constraints
    const auto& c = m_machine.getConstraints();
    std::cout << "\n  Machine Constraints:\n";
    std::cout << "    Gantry Speed:  [" << c.gantryRotationSpeed[0] << ", " << c.gantryRotationSpeed[1] << "] deg/s\n";
    std::cout << "    Leaf Speed:    [" << c.leafSpeed[0] << ", " << c.leafSpeed[1] << "] mm/s\n";
    std::cout << "    MU Rate:       [" << c.monitorUnitRate[0] << ", " << c.monitorUnitRate[1] << "] MU/s\n";

    std::cout << "==============================\n";
}

} // namespace optirad
