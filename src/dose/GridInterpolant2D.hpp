#pragma once

#include <vector>
#include <cmath>
#include <algorithm>

namespace optirad {

/**
 * Simple 2D bilinear interpolation on a regular grid.
 * Used to evaluate convolved kernel profiles at arbitrary lateral positions.
 * Analogous to MATLAB's griddedInterpolant.
 */
class GridInterpolant2D {
public:
    GridInterpolant2D() = default;

    /**
     * Construct interpolant from regular grid data.
     * @param xMin, xMax  Range of x coordinates
     * @param yMin, yMax  Range of y coordinates
     * @param nx, ny      Number of grid points in x and y
     * @param data        Row-major data array of size nx * ny
     */
    GridInterpolant2D(double xMin, double xMax, size_t nx,
                      double yMin, double yMax, size_t ny,
                      const std::vector<double>& data)
        : m_xMin(xMin), m_xMax(xMax), m_nx(nx)
        , m_yMin(yMin), m_yMax(yMax), m_ny(ny)
        , m_data(data)
    {
        m_dx = (nx > 1) ? (xMax - xMin) / static_cast<double>(nx - 1) : 1.0;
        m_dy = (ny > 1) ? (yMax - yMin) / static_cast<double>(ny - 1) : 1.0;
    }

    /// Set grid data after construction
    void setGrid(double xMin, double xMax, size_t nx,
                 double yMin, double yMax, size_t ny,
                 const std::vector<double>& data) {
        m_xMin = xMin; m_xMax = xMax; m_nx = nx;
        m_yMin = yMin; m_yMax = yMax; m_ny = ny;
        m_data = data;
        m_dx = (nx > 1) ? (xMax - xMin) / static_cast<double>(nx - 1) : 1.0;
        m_dy = (ny > 1) ? (yMax - yMin) / static_cast<double>(ny - 1) : 1.0;
    }

    /// Evaluate at a single point using bilinear interpolation.
    /// Returns 0 for out-of-bounds queries (extrapolation = 0).
    double operator()(double x, double y) const {
        // Normalize to grid coordinates
        double fx = (x - m_xMin) / m_dx;
        double fy = (y - m_yMin) / m_dy;

        // Check bounds - return 0 outside
        if (fx < 0.0 || fx > static_cast<double>(m_nx - 1) ||
            fy < 0.0 || fy > static_cast<double>(m_ny - 1)) {
            return 0.0;
        }

        size_t ix = static_cast<size_t>(fx);
        size_t iy = static_cast<size_t>(fy);

        // Clamp to valid range
        if (ix >= m_nx - 1) ix = m_nx - 2;
        if (iy >= m_ny - 1) iy = m_ny - 2;

        double tx = fx - static_cast<double>(ix);
        double ty = fy - static_cast<double>(iy);

        // Bilinear interpolation
        double v00 = m_data[ix * m_ny + iy];
        double v10 = m_data[(ix + 1) * m_ny + iy];
        double v01 = m_data[ix * m_ny + (iy + 1)];
        double v11 = m_data[(ix + 1) * m_ny + (iy + 1)];

        return (1.0 - tx) * (1.0 - ty) * v00 +
               tx         * (1.0 - ty) * v10 +
               (1.0 - tx) * ty         * v01 +
               tx         * ty         * v11;
    }

    /// Evaluate at multiple points
    std::vector<double> evaluate(const std::vector<double>& xs, const std::vector<double>& ys) const {
        size_t n = xs.size();
        std::vector<double> result(n);
        for (size_t i = 0; i < n; ++i) {
            result[i] = (*this)(xs[i], ys[i]);
        }
        return result;
    }

    bool isValid() const { return !m_data.empty() && m_nx > 0 && m_ny > 0; }

private:
    double m_xMin = 0.0, m_xMax = 1.0;
    double m_yMin = 0.0, m_yMax = 1.0;
    size_t m_nx = 0, m_ny = 0;
    double m_dx = 1.0, m_dy = 1.0;
    std::vector<double> m_data;
};

} // namespace optirad
