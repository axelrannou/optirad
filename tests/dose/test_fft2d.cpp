#include <gtest/gtest.h>
#include "dose/FFT2D.hpp"
#include <cmath>
#include <numeric>

namespace optirad::tests {

// ============================================================================
// nextPow2
// ============================================================================

TEST(FFT2DTest, NextPow2) {
    EXPECT_EQ(FFT2D::nextPow2(1), 1u);
    EXPECT_EQ(FFT2D::nextPow2(2), 2u);
    EXPECT_EQ(FFT2D::nextPow2(3), 4u);
    EXPECT_EQ(FFT2D::nextPow2(4), 4u);
    EXPECT_EQ(FFT2D::nextPow2(5), 8u);
    EXPECT_EQ(FFT2D::nextPow2(1023), 1024u);
    EXPECT_EQ(FFT2D::nextPow2(1024), 1024u);
}

// ============================================================================
// 2D FFT forward → inverse should reconstruct the original signal
// ============================================================================

TEST(FFT2DTest, ForwardInverseRoundTrip) {
    size_t n = 8;
    std::vector<FFT2D::Complex> data(n * n);
    for (size_t i = 0; i < n * n; ++i)
        data[i] = {static_cast<double>(i + 1), 0.0};

    // Save original
    auto original = data;

    // Forward FFT
    FFT2D::fft2(data, n, n, false);

    // Inverse FFT
    FFT2D::fft2(data, n, n, true);

    // Should match original
    for (size_t i = 0; i < n * n; ++i) {
        EXPECT_NEAR(data[i].real(), original[i].real(), 1e-10)
            << "Real mismatch at index " << i;
        EXPECT_NEAR(data[i].imag(), original[i].imag(), 1e-10)
            << "Imag mismatch at index " << i;
    }
}

// ============================================================================
// 2D convolution with Dirac delta → identity
// ============================================================================

TEST(FFT2DTest, ConvolveWithDelta) {
    // Convolving any signal with a centered delta should return the signal
    size_t n = 5;
    std::vector<double> signal(n * n);
    for (size_t i = 0; i < n * n; ++i)
        signal[i] = static_cast<double>(i + 1);

    // 1x1 delta
    std::vector<double> delta = {1.0};
    auto result = FFT2D::convolve2D(signal, n, n, delta, 1, 1);

    // Result should be the same as input (full conv with 1x1 = same size)
    ASSERT_EQ(result.size(), n * n);
    for (size_t i = 0; i < n * n; ++i) {
        EXPECT_NEAR(result[i], signal[i], 1e-10) << "at index " << i;
    }
}

// ============================================================================
// 2D convolution commutativity: a*b == b*a
// ============================================================================

TEST(FFT2DTest, ConvolutionCommutative) {
    std::vector<double> a = {1, 2, 3, 4};     // 2x2
    std::vector<double> b = {5, 6, 7, 8, 9};  // 1x5

    auto ab = FFT2D::convolve2D(a, 2, 2, b, 1, 5);
    auto ba = FFT2D::convolve2D(b, 1, 5, a, 2, 2);

    // Both produce 2x6 = 12 values
    ASSERT_EQ(ab.size(), ba.size());
    for (size_t i = 0; i < ab.size(); ++i) {
        EXPECT_NEAR(ab[i], ba[i], 1e-10) << "at index " << i;
    }
}

// ============================================================================
// convolve2DSame returns same size as first input
// ============================================================================

TEST(FFT2DTest, ConvoleSameSize) {
    size_t aR = 7, aC = 9;
    size_t bR = 3, bC = 3;

    std::vector<double> a(aR * aC, 1.0);
    std::vector<double> b(bR * bC, 1.0);

    auto result = FFT2D::convolve2DSame(a, aR, aC, b, bR, bC);
    EXPECT_EQ(result.size(), aR * aC);
}

// ============================================================================
// Known 1D convolution via FFT2D (use 1xN grids)
// ============================================================================

TEST(FFT2DTest, Known1DConvolution) {
    // [1, 2, 3] * [1, 1] = [1, 3, 5, 3]
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {1.0, 1.0};
    auto result = FFT2D::convolve2D(a, 1, 3, b, 1, 2);

    ASSERT_EQ(result.size(), 4u); // 3+2-1 = 4 (1x4)
    EXPECT_NEAR(result[0], 1.0, 1e-10);
    EXPECT_NEAR(result[1], 3.0, 1e-10);
    EXPECT_NEAR(result[2], 5.0, 1e-10);
    EXPECT_NEAR(result[3], 3.0, 1e-10);
}

} // namespace optirad::tests
