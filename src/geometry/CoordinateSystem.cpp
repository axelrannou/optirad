#include "CoordinateSystem.hpp"

namespace optirad {

Vec3 CoordinateSystem::patientToWorld(const Vec3& patient) const {
    return patient; // TODO: Apply transformation matrix
}

Vec3 CoordinateSystem::worldToPatient(const Vec3& world) const {
    return world; // TODO: Apply inverse transformation
}

Vec3 CoordinateSystem::voxelToWorld(int i, int j, int k) const {
    return {
        m_origin[0] + i * m_spacing[0],
        m_origin[1] + j * m_spacing[1],
        m_origin[2] + k * m_spacing[2]
    };
}

void CoordinateSystem::worldToVoxel(const Vec3& world, int& i, int& j, int& k) const {
    // Check for zero spacing to prevent division by zero
    if (std::abs(m_spacing[0]) < 1e-10 || std::abs(m_spacing[1]) < 1e-10 || std::abs(m_spacing[2]) < 1e-10) {
        i = j = k = 0;
        return;
    }
    i = static_cast<int>((world[0] - m_origin[0]) / m_spacing[0]);
    j = static_cast<int>((world[1] - m_origin[1]) / m_spacing[1]);
    k = static_cast<int>((world[2] - m_origin[2]) / m_spacing[2]);
}

void CoordinateSystem::setOrigin(const Vec3& origin) { m_origin = origin; }
void CoordinateSystem::setSpacing(const Vec3& spacing) { m_spacing = spacing; }

} // namespace optirad
