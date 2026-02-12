#include "StructureSet.hpp"
#include "Grid.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>

namespace optirad {

void StructureSet::addStructure(std::unique_ptr<Structure> structure) {
    m_structures.push_back(std::move(structure));
}

const Structure* StructureSet::getStructure(size_t index) const {
    if (index >= m_structures.size()) return nullptr;
    return m_structures[index].get();
}

Structure* StructureSet::getStructure(size_t index) {
    if (index >= m_structures.size()) return nullptr;
    return m_structures[index].get();
}

const Structure* StructureSet::getStructureByName(const std::string& name) const {
    for (const auto& s : m_structures) {
        if (s->getName() == name) return s.get();
    }
    return nullptr;
}

void StructureSet::rasterizeContours(const Grid& ctGrid) {
    Logger::info("Rasterizing structure contours...");
    
    for (auto& structure : m_structures) {
        if (!structure) continue;
        
        structure->rasterizeContours(ctGrid);
        
        Logger::info("  - Rasterizing: " + structure->getName() + 
                     " -> " + std::to_string(structure->getVoxelIndices().size()) + " voxels");
    }
    
    Logger::info("Rasterization complete.");
}

} // namespace optirad
