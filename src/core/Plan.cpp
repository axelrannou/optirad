#include "Plan.hpp"
#include "PatientData.hpp"
#include "geometry/Volume.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "utils/Logger.hpp"

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

    // Collect all voxel indices from PTV structures only
    std::set<size_t> ptvVoxels;
    int ptvCount = 0;
    
    for (size_t i = 0; i < structureSet->getCount(); ++i) {
        const auto* structure = structureSet->getStructure(i);
        // I will fix it later but to compare with matRad we keep all "targets" type
        if (structure && structure->getType() == "PTV" || structure->getType() == "GTV" || structure->getType() == "CTV") {
            ptvCount++;
            // Get voxel indices for this PTV structure
            const auto& voxels = structure->getVoxelIndices();
            if (!voxels.empty()) {
                ptvVoxels.insert(voxels.begin(), voxels.end());
            } else {
                Logger::warn("PTV structure '" + structure->getName() + 
                           "' has no voxel indices (contours not rasterized yet)");
            }
        }
    }

    // Error out if no PTV structures found
    if (ptvCount == 0) {
        Logger::error("No PTV structures found in structure set");
        Logger::error("Available structures:");
        for (size_t i = 0; i < structureSet->getCount(); ++i) {
            const auto* s = structureSet->getStructure(i);
            if (s) {
                Logger::error("  - " + s->getName() + " (type: " + s->getType() + ")");
            }
        }
        Logger::error("Cannot compute isocenter without PTV structures");
        return {0.0, 0.0, 0.0};
    }

    // Error out if PTV structures have no voxel indices
    if (ptvVoxels.empty()) {
        Logger::error("PTV structures found but have no voxel indices");
        Logger::error("You must rasterize structure contours before computing isocenter");
        return {0.0, 0.0, 0.0};
    }

    // Convert voxel indices to world coordinates and compute center of gravity
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    
    size_t totalVoxels = dims[0] * dims[1] * dims[2];
    for (size_t voxelIdx : ptvVoxels) {
        if (voxelIdx >= totalVoxels) {
            Logger::warn("Plan::computeIsoCenter: voxel index " + std::to_string(voxelIdx) + 
                        " out of bounds (valid range: 0 to " + std::to_string(totalVoxels - 1) + "), skipping");
            continue;
        }
        
        // Convert linear index to 3D coordinates (column-major: dims = [ny, nx, nz])
        // index = i + j*ny + k*ny*nx where i=y-index, j=x-index, k=z-index
        size_t k = voxelIdx / (dims[0] * dims[1]);  // z-index (slice)
        size_t j = (voxelIdx % (dims[0] * dims[1])) / dims[0];  // x-index (column)
        size_t i = voxelIdx % dims[0];  // y-index (row)
        
        // Convert to world coordinates using the Grid's transformation
        Vec3 worldPos = grid.voxelToPatient({static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
        
        sumX += worldPos[0];
        sumY += worldPos[1];
        sumZ += worldPos[2];
    }

    size_t numVoxels = ptvVoxels.size();
    std::array<double, 3> isoCenter;
    isoCenter[0] = sumX / numVoxels;
    isoCenter[1] = sumY / numVoxels;
    isoCenter[2] = sumZ / numVoxels;

    Logger::info("Computed isocenter from " + std::to_string(numVoxels) + " PTV voxels across " + 
                std::to_string(ptvCount) + " PTV structure(s)");

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
