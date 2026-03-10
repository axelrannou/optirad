#include "ObjectiveFunction.hpp"

namespace optirad {

void ObjectiveFunction::setWeight(double weight) { m_weight = weight; }
double ObjectiveFunction::getWeight() const { return m_weight; }
void ObjectiveFunction::setStructure(const Structure* structure) { m_structure = structure; }

void ObjectiveFunction::setVoxelIndices(const std::vector<size_t>& indices) {
    m_mappedIndices = indices;
}

const std::vector<size_t>& ObjectiveFunction::getActiveIndices() const {
    if (!m_mappedIndices.empty()) return m_mappedIndices;
    // Fallback to structure indices (WARNING: may be CT-grid indexed!)
    static const std::vector<size_t> empty;
    if (!m_structure) return empty;
    return m_structure->getVoxelIndices();
}

} // namespace optirad
