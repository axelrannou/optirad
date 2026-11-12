#include "ObjectiveFunction.hpp"
#include <sstream>

namespace optirad {

std::string ObjectiveFunction::describe() const {
    std::ostringstream out;
    out << (m_structure ? m_structure->getName() : "<unbound>") << " / " << getName();
    DoseObjectiveSpec spec;
    if (gpuDescription(spec)) {
        out << " dose=" << spec.threshold << " Gy";
        if (spec.type == 3 || spec.type == 4) out << " volume=" << 100.0 * spec.volumeFraction << "%";
    }
    out << " weight=" << m_weight << " voxels=" << getActiveIndices().size();
    return out.str();
}

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
