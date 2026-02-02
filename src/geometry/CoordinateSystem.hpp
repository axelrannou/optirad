#pragma once

#include "MathUtils.hpp"

namespace optirad {

class CoordinateSystem {
public:
    // Transform from DICOM patient coordinates to internal coordinates
    Vec3 patientToWorld(const Vec3& patient) const;
    Vec3 worldToPatient(const Vec3& world) const;

    // Transform from voxel indices to world coordinates
    Vec3 voxelToWorld(int i, int j, int k) const;
    void worldToVoxel(const Vec3& world, int& i, int& j, int& k) const;

    void setOrigin(const Vec3& origin);
    void setSpacing(const Vec3& spacing);

private:
    Vec3 m_origin = {0.0, 0.0, 0.0};
    Vec3 m_spacing = {1.0, 1.0, 1.0};
};

} // namespace optirad
