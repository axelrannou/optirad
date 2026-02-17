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
    // Protect against division by zero - return zero vector if norm is too small
    if (std::abs(n) < 1e-14) {
        return {0.0, 0.0, 0.0};
    }
    return {v[0]/n, v[1]/n, v[2]/n};
}

// Matrix transpose for 3x3
inline Mat3 transpose(const Mat3& m) {
    return Mat3(
        m.m[0][0], m.m[1][0], m.m[2][0],
        m.m[0][1], m.m[1][1], m.m[2][1],
        m.m[0][2], m.m[1][2], m.m[2][2]
    );
}

// Matrix-matrix multiplication
inline Mat3 matMul(const Mat3& a, const Mat3& b) {
    Mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result.m[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                result.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return result;
}

/**
 * Compute the combined gantry + couch rotation matrix.
 * 
 * Gantry rotation: active counter-clockwise rotation around the Z-axis (LPS superior).
 * Couch rotation: active counter-clockwise rotation around the Y-axis (LPS anterior-posterior).
 * Combined: R = R_Couch * R_Gantry
 * 
 * Angles are in degrees.
 */
inline Mat3 getRotationMatrix(double gantryAngleDeg, double couchAngleDeg) {
    constexpr double deg2rad = M_PI / 180.0;
    double gantryRad = gantryAngleDeg * deg2rad;
    double couchRad  = couchAngleDeg * deg2rad;

    double cg = std::cos(gantryRad);
    double sg = std::sin(gantryRad);
    double cc = std::cos(couchRad);
    double sc = std::sin(couchRad);

    // Gantry rotation around Z-axis
    Mat3 Rg(
         cg, -sg, 0.0,
         sg,  cg, 0.0,
        0.0, 0.0, 1.0
    );

    // Couch rotation around Y-axis
    Mat3 Rc(
         cc, 0.0, sc,
        0.0, 1.0, 0.0,
        -sc, 0.0, cc
    );

    // Combined: R = R_Couch * R_Gantry
    return matMul(Rc, Rg);
}

// Vector addition/subtraction
inline Vec3 vecAdd(const Vec3& a, const Vec3& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

inline Vec3 vecSub(const Vec3& a, const Vec3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

inline Vec3 vecScale(const Vec3& v, double s) {
    return {v[0]*s, v[1]*s, v[2]*s};
}

// Multiply row vector by matrix transpose: v * M^T  (equivalent to (M * v^T)^T for column vectors)
// This matches matRad's convention: coords * rotMat_vectors_T
inline Vec3 vecMulMatTranspose(const Vec3& v, const Mat3& m) {
    return {
        v[0]*m.m[0][0] + v[1]*m.m[0][1] + v[2]*m.m[0][2],
        v[0]*m.m[1][0] + v[1]*m.m[1][1] + v[2]*m.m[1][2],
        v[0]*m.m[2][0] + v[1]*m.m[2][1] + v[2]*m.m[2][2]
    };
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
