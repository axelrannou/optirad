#pragma once

#include "MathUtils.hpp"
#include <array>

namespace optirad {

class Grid {
public:
    void setDimensions(size_t nx, size_t ny, size_t nz);
    void setSpacing(double dx, double dy, double dz);
    void setOrigin(const Vec3& origin);

    std::array<size_t, 3> getDimensions() const;
    Vec3 getSpacing() const;
    Vec3 getOrigin() const;
    size_t getNumVoxels() const;

private:
    std::array<size_t, 3> m_dimensions = {0, 0, 0};
    Vec3 m_spacing = {1.0, 1.0, 1.0};
    Vec3 m_origin = {0.0, 0.0, 0.0};
};

} // namespace optirad
