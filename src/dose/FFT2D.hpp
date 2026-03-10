#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace optirad {

/**
 * Lightweight 2D FFT for convolution in photon pencil beam dose calculation.
 * Uses Cooley-Tukey radix-2 algorithm. Input sizes are automatically zero-padded to next power of 2.
 */
class FFT2D {
public:
    using Complex = std::complex<double>;

    /// Perform 2D FFT in-place on data of size rows x cols (must be power of 2)
    static void fft2(std::vector<Complex>& data, size_t rows, size_t cols, bool inverse = false) {
        // FFT along rows
        for (size_t r = 0; r < rows; ++r) {
            std::vector<Complex> row(cols);
            for (size_t c = 0; c < cols; ++c)
                row[c] = data[r * cols + c];
            fft1d(row, inverse);
            for (size_t c = 0; c < cols; ++c)
                data[r * cols + c] = row[c];
        }
        // FFT along columns
        for (size_t c = 0; c < cols; ++c) {
            std::vector<Complex> col(rows);
            for (size_t r = 0; r < rows; ++r)
                col[r] = data[r * cols + c];
            fft1d(col, inverse);
            for (size_t r = 0; r < rows; ++r)
                data[r * cols + c] = col[r];
        }
    }

    /// Perform 2D convolution of two real-valued grids using FFT
    /// a and b are row-major, size rows x cols. Result is same size.
    static std::vector<double> convolve2D(
        const std::vector<double>& a, size_t aRows, size_t aCols,
        const std::vector<double>& b, size_t bRows, size_t bCols)
    {
        // Output size for full convolution
        size_t outRows = aRows + bRows - 1;
        size_t outCols = aCols + bCols - 1;

        // Pad to next power of 2
        size_t fftRows = nextPow2(outRows);
        size_t fftCols = nextPow2(outCols);

        // Zero-pad and convert to complex
        std::vector<Complex> A(fftRows * fftCols, {0.0, 0.0});
        std::vector<Complex> B(fftRows * fftCols, {0.0, 0.0});

        for (size_t r = 0; r < aRows; ++r)
            for (size_t c = 0; c < aCols; ++c)
                A[r * fftCols + c] = {a[r * aCols + c], 0.0};

        for (size_t r = 0; r < bRows; ++r)
            for (size_t c = 0; c < bCols; ++c)
                B[r * fftCols + c] = {b[r * bCols + c], 0.0};

        // Forward FFT
        fft2(A, fftRows, fftCols, false);
        fft2(B, fftRows, fftCols, false);

        // Point-wise multiply
        for (size_t i = 0; i < A.size(); ++i)
            A[i] *= B[i];

        // Inverse FFT
        fft2(A, fftRows, fftCols, true);

        // Extract real part at original output size
        std::vector<double> result(outRows * outCols);
        for (size_t r = 0; r < outRows; ++r)
            for (size_t c = 0; c < outCols; ++c)
                result[r * outCols + c] = A[r * fftCols + c].real();

        return result;
    }

    /// Convolve two grids and return result cropped to 'same' size (centered on a)
    static std::vector<double> convolve2DSame(
        const std::vector<double>& a, size_t aRows, size_t aCols,
        const std::vector<double>& b, size_t bRows, size_t bCols)
    {
        auto full = convolve2D(a, aRows, aCols, b, bRows, bCols);
        size_t fullRows = aRows + bRows - 1;
        size_t fullCols = aCols + bCols - 1;

        // Crop to center (same size as a)
        size_t offsetR = bRows / 2;
        size_t offsetC = bCols / 2;

        std::vector<double> result(aRows * aCols);
        for (size_t r = 0; r < aRows; ++r)
            for (size_t c = 0; c < aCols; ++c)
                result[r * aCols + c] = full[(r + offsetR) * fullCols + (c + offsetC)];

        return result;
    }

    /// Convolve two grids using FFT with explicit zero-padding to padSize×padSize.
    /// Replicates MATLAB: real(ifft2(fft2(a,N,N).*fft2(b,N,N))) where N=padSize.
    /// padSize must be >= aRows + bRows - 1 for correct linear convolution.
    static std::vector<double> convolve2DPadded(
        const std::vector<double>& a, size_t aRows, size_t aCols,
        const std::vector<double>& b, size_t bRows, size_t bCols,
        size_t padSize)
    {
        // Pad to power of 2 for Cooley-Tukey FFT
        size_t fftN = nextPow2(padSize);

        std::vector<Complex> A(fftN * fftN, {0.0, 0.0});
        std::vector<Complex> B(fftN * fftN, {0.0, 0.0});

        // Zero-pad a into top-left corner
        for (size_t r = 0; r < aRows; ++r)
            for (size_t c = 0; c < aCols; ++c)
                A[r * fftN + c] = {a[r * aCols + c], 0.0};

        // Zero-pad b into top-left corner
        for (size_t r = 0; r < bRows; ++r)
            for (size_t c = 0; c < bCols; ++c)
                B[r * fftN + c] = {b[r * bCols + c], 0.0};

        // Forward FFT
        fft2(A, fftN, fftN, false);
        fft2(B, fftN, fftN, false);

        // Point-wise multiply
        for (size_t i = 0; i < A.size(); ++i)
            A[i] *= B[i];

        // Inverse FFT
        fft2(A, fftN, fftN, true);

        // Extract padSize × padSize result
        std::vector<double> result(padSize * padSize);
        for (size_t r = 0; r < padSize; ++r)
            for (size_t c = 0; c < padSize; ++c)
                result[r * padSize + c] = A[r * fftN + c].real();

        return result;
    }

    static size_t nextPow2(size_t n) {
        size_t p = 1;
        while (p < n) p <<= 1;
        return p;
    }

private:
    /// In-place 1D FFT (Cooley-Tukey radix-2 DIT). data.size() must be power of 2.
    static void fft1d(std::vector<Complex>& data, bool inverse) {
        size_t n = data.size();
        if (n <= 1) return;

        // Bit-reversal permutation
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j) std::swap(data[i], data[j]);
        }

        // Butterfly operations
        for (size_t len = 2; len <= n; len <<= 1) {
            double angle = 2.0 * M_PI / static_cast<double>(len) * (inverse ? -1.0 : 1.0);
            Complex wlen(std::cos(angle), std::sin(angle));

            for (size_t i = 0; i < n; i += len) {
                Complex w(1.0, 0.0);
                for (size_t j = 0; j < len / 2; ++j) {
                    Complex u = data[i + j];
                    Complex v = w * data[i + j + len / 2];
                    data[i + j] = u + v;
                    data[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        // Scale for inverse
        if (inverse) {
            double invN = 1.0 / static_cast<double>(n);
            for (auto& x : data)
                x *= invN;
        }
    }
};

} // namespace optirad
