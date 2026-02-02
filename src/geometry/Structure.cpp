#include "Structure.hpp"

namespace optirad {

void Structure::setName(const std::string& name) { m_name = name; }
void Structure::setColor(float r, float g, float b) {
    m_color[0] = r; m_color[1] = g; m_color[2] = b;
}
void Structure::setType(const std::string& type) { m_type = type; }

const std::string& Structure::getName() const { return m_name; }
const std::string& Structure::getType() const { return m_type; }

void Structure::setMask(const std::vector<bool>& mask) { m_mask = mask; }
const std::vector<bool>& Structure::getMask() const { return m_mask; }

std::vector<size_t> Structure::getVoxelIndices() const {
    std::vector<size_t> indices;
    for (size_t i = 0; i < m_mask.size(); ++i) {
        if (m_mask[i]) indices.push_back(i);
    }
    return indices;
}

} // namespace optirad
