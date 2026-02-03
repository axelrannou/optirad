#pragma once

#include <array>

namespace optirad {

using Vec3 = std::array<double, 3>;

namespace MathUtils {
    double dot(const Vec3& a, const Vec3& b);
    Vec3 cross(const Vec3& a, const Vec3& b);
    double norm(const Vec3& v);
    Vec3 normalize(const Vec3& v);
}

} // namespace optirad
