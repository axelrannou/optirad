#pragma once

#include <array>
#include <cmath>

namespace optirad {

// 3D vector using std::array
using Vec3 = std::array<double, 3>;

// 3x3 matrix for rotations and transforms
struct Mat3 {
    double m[3][3];
    
    Mat3() : m{{1,0,0},{0,1,0},{0,0,1}} {}
    
    Mat3(double m00, double m01, double m02,
         double m10, double m11, double m12,
         double m20, double m21, double m22)
        : m{{m00, m01, m02}, {m10, m11, m12}, {m20, m21, m22}} {}
    
    Vec3 operator*(const Vec3& v) const {
        return {
            m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
            m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
            m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
        };
    }
};

// Vector operations - all inline for performance
inline double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    if (std::abs(n) < 1e-14) {
        return {0.0, 0.0, 0.0};  // Return zero vector for zero-length input
    }
    return {v[0]/n, v[1]/n, v[2]/n};
}

// Matrix inverse for 3x3
inline Mat3 inverse(const Mat3& m) {
    double det = m.m[0][0]*(m.m[1][1]*m.m[2][2] - m.m[1][2]*m.m[2][1])
               - m.m[0][1]*(m.m[1][0]*m.m[2][2] - m.m[1][2]*m.m[2][0])
               + m.m[0][2]*(m.m[1][0]*m.m[2][1] - m.m[1][1]*m.m[2][0]);
    
    // Check for singular matrix
    if (std::abs(det) < 1e-14) {
        // Return identity matrix for singular matrix
        return Mat3();
    }
    
    double invDet = 1.0 / det;
    
    Mat3 result;
    result.m[0][0] = (m.m[1][1]*m.m[2][2] - m.m[1][2]*m.m[2][1]) * invDet;
    result.m[0][1] = (m.m[0][2]*m.m[2][1] - m.m[0][1]*m.m[2][2]) * invDet;
    result.m[0][2] = (m.m[0][1]*m.m[1][2] - m.m[0][2]*m.m[1][1]) * invDet;
    result.m[1][0] = (m.m[1][2]*m.m[2][0] - m.m[1][0]*m.m[2][2]) * invDet;
    result.m[1][1] = (m.m[0][0]*m.m[2][2] - m.m[0][2]*m.m[2][0]) * invDet;
    result.m[1][2] = (m.m[0][2]*m.m[1][0] - m.m[0][0]*m.m[1][2]) * invDet;
    result.m[2][0] = (m.m[1][0]*m.m[2][1] - m.m[1][1]*m.m[2][0]) * invDet;
    result.m[2][1] = (m.m[0][1]*m.m[2][0] - m.m[0][0]*m.m[2][1]) * invDet;
    result.m[2][2] = (m.m[0][0]*m.m[1][1] - m.m[0][1]*m.m[1][0]) * invDet;
    return result;
}

} // namespace optirad
