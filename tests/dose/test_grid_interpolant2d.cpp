#include <gtest/gtest.h>
#include "dose/GridInterpolant2D.hpp"
#include <cmath>

namespace optirad::tests {

// ============================================================================
// Basic construction and grid-point queries
// ============================================================================

TEST(GridInterpolant2DTest, EmptyGridReturnsZero) {
    // Grid with 0 data evaluates to 0
    GridInterpolant2D interp(0, 1, 2, 0, 1, 2, {0.0, 0.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(interp(0.5, 0.5), 0.0);
}

TEST(GridInterpolant2DTest, ConstantGrid) {
    // 3x3 grid all 5.0 from [-1,1] x [-1,1]
    std::vector<double> data(9, 5.0);
    GridInterpolant2D interp(-1, 1, 3, -1, 1, 3, data);

    EXPECT_DOUBLE_EQ(interp(0.0, 0.0), 5.0);
    EXPECT_DOUBLE_EQ(interp(-1.0, -1.0), 5.0);
    EXPECT_DOUBLE_EQ(interp(1.0, 1.0), 5.0);
    EXPECT_DOUBLE_EQ(interp(0.5, 0.5), 5.0);
}

TEST(GridInterpolant2DTest, LinearInterpolationX) {
    // 2x2 grid: x ∈ [0, 1], y ∈ [0, 1]
    // data constant along y: row0={0, 0}, row1={10, 10}
    // At x=0.5, y=0.5 => bilinear gives 5.0
    std::vector<double> data = {0.0, 0.0, 10.0, 10.0};
    GridInterpolant2D interp(0, 1, 2, 0, 1, 2, data);

    EXPECT_NEAR(interp(0.0, 0.5), 0.0, 1e-10);
    EXPECT_NEAR(interp(1.0, 0.5), 10.0, 1e-10);
    EXPECT_NEAR(interp(0.5, 0.5), 5.0, 1e-10);
    EXPECT_NEAR(interp(0.25, 0.5), 2.5, 1e-10);
}

TEST(GridInterpolant2DTest, BilinearInterpolation) {
    // 2x2 grid: x ∈ [0, 1], y ∈ [0, 1]
    // data[0*2+0]=0 data[0*2+1]=1 data[1*2+0]=2 data[1*2+1]=3
    // At (0.5, 0.5): bilinear = (0+1+2+3)/4 = 1.5
    std::vector<double> data = {0.0, 1.0, 2.0, 3.0};
    GridInterpolant2D interp(0, 1, 2, 0, 1, 2, data);

    EXPECT_NEAR(interp(0.5, 0.5), 1.5, 1e-10);
    EXPECT_DOUBLE_EQ(interp(0.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(interp(1.0, 1.0), 3.0);
}

// ============================================================================
// Out-of-bounds returns 0
// ============================================================================

TEST(GridInterpolant2DTest, OutOfBoundsReturnsZero) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    GridInterpolant2D interp(0, 1, 2, 0, 1, 2, data);

    EXPECT_DOUBLE_EQ(interp(-0.1, 0.5), 0.0);
    EXPECT_DOUBLE_EQ(interp(1.1, 0.5), 0.0);
    EXPECT_DOUBLE_EQ(interp(0.5, -0.1), 0.0);
    EXPECT_DOUBLE_EQ(interp(0.5, 1.1), 0.0);
}

// ============================================================================
// setGrid after construction
// ============================================================================

TEST(GridInterpolant2DTest, SetGridOverrides) {
    GridInterpolant2D interp;
    std::vector<double> data = {10.0, 20.0, 30.0, 40.0};
    interp.setGrid(0, 1, 2, 0, 1, 2, data);

    EXPECT_DOUBLE_EQ(interp(0.0, 0.0), 10.0);
    EXPECT_DOUBLE_EQ(interp(1.0, 1.0), 40.0);
}

// ============================================================================
// Larger grid symmetry check
// ============================================================================

TEST(GridInterpolant2DTest, SymmetricGaussian) {
    // Create a symmetric Gaussian pattern on a 5x5 grid from [-2,2] x [-2,2]
    size_t n = 5;
    std::vector<double> data(n * n);
    double sigma = 1.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double x = -2.0 + i;
            double y = -2.0 + j;
            data[i * n + j] = std::exp(-(x * x + y * y) / (2 * sigma * sigma));
        }
    }

    GridInterpolant2D interp(-2, 2, 5, -2, 2, 5, data);

    // Gaussian is symmetric: f(1,0) == f(-1,0) == f(0,1)
    double v_pos = interp(1.0, 0.0);
    double v_neg = interp(-1.0, 0.0);
    double v_perp = interp(0.0, 1.0);
    EXPECT_NEAR(v_pos, v_neg, 1e-10);
    EXPECT_NEAR(v_pos, v_perp, 1e-10);

    // Center should be maximum
    double v_center = interp(0.0, 0.0);
    EXPECT_GT(v_center, v_pos);
}

} // namespace optirad::tests
