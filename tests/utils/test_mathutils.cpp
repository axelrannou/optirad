#include <gtest/gtest.h>
#include "geometry/MathUtils.hpp"
#include <cmath>

namespace optirad::tests {

class MathUtilsTest : public ::testing::Test {
protected:
    // Test mathematical utility functions
};

TEST(MathUtilsTest, NormalizeValidVector) {
    Vec3 v{3.0, 4.0, 0.0};  // 3-4-5 triangle
    auto normalized = normalize(v);
    
    double length = norm(normalized);
    EXPECT_NEAR(length, 1.0, 1e-10);
}

TEST(MathUtilsTest, NormalizeZeroVector) {
    Vec3 zero{0.0, 0.0, 0.0};
    auto normalized = normalize(zero);
    // Should return zero vector for zero-length input
    EXPECT_EQ(normalized[0], 0.0);
    EXPECT_EQ(normalized[1], 0.0);
    EXPECT_EQ(normalized[2], 0.0);
}

TEST(MathUtilsTest, DotProductOrthogonal) {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    
    double dotProd = dot(a, b);
    EXPECT_DOUBLE_EQ(dotProd, 0.0);
}

TEST(MathUtilsTest, DotProductParallel) {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{2.0, 0.0, 0.0};
    
    double dotProd = dot(a, b);
    EXPECT_DOUBLE_EQ(dotProd, 2.0);
}

TEST(MathUtilsTest, CrossProductOrthogonal) {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    
    auto crossProd = cross(a, b);
    // Should point in z direction
    EXPECT_DOUBLE_EQ(crossProd[2], 1.0);
    EXPECT_DOUBLE_EQ(crossProd[0], 0.0);
    EXPECT_DOUBLE_EQ(crossProd[1], 0.0);
}

TEST(MathUtilsTest, VectorNorm) {
    Vec3 v{3.0, 4.0, 0.0};
    double len = norm(v);
    EXPECT_DOUBLE_EQ(len, 5.0);  // 3-4-5 triangle
}

TEST(MathUtilsTest, Mat3MatrixMultiplication) {
    // Identity matrix
    Mat3 identity;
    Vec3 v{1.0, 2.0, 3.0};
    
    Vec3 result = identity * v;
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
}

TEST(MathUtilsTest, Mat3InverseOfIdentity) {
    Mat3 identity;
    Mat3 inv = inverse(identity);
    
    // Inverse of identity should be identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(inv.m[i][j], expected, 1e-10);
        }
    }
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
