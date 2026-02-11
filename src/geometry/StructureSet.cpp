#include "StructureSet.hpp"
#include "Grid.hpp"
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
    std::cout << "Rasterizing structure contours...\n";
    
    for (auto& structure : m_structures) {
        if (!structure) continue;
        
        std::cout << "  - Rasterizing: " << structure->getName() << " ... ";
        
        structure->rasterizeContours(ctGrid);
        
        std::cout << structure->getVoxelIndices().size() << " voxels\n";
    }
    
    std::cout << "Rasterization complete.\n";
}

} // namespace optirad
