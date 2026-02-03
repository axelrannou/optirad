#include "StructureSet.hpp"

namespace optirad {

void StructureSet::addStructure(const Structure& structure) {
    m_structures.push_back(structure);
}

size_t StructureSet::getNumStructures() const { return m_structures.size(); }

const Structure& StructureSet::getStructure(size_t index) const {
    return m_structures[index];
}

std::optional<const Structure*> StructureSet::findByName(const std::string& name) const {
    for (const auto& s : m_structures) {
        if (s.getName() == name) return &s;
    }
    return std::nullopt;
}

} // namespace optirad
