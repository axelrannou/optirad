#include "DoseMatrix.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace optirad {

void DoseMatrix::setGrid(const Grid& grid) { m_grid = grid; }

void DoseMatrix::allocate() {
    m_data.resize(m_grid.getNumVoxels(), 0.0);
}

double& DoseMatrix::at(size_t i, size_t j, size_t k) {
    auto dims = m_grid.getDimensions();
    if (i >= dims[0] || j >= dims[1] || k >= dims[2]) {
        throw std::out_of_range("DoseMatrix::at: index out of bounds");
    }
    return m_data[i + j * dims[0] + k * dims[0] * dims[1]];
}

const double& DoseMatrix::at(size_t i, size_t j, size_t k) const {
    auto dims = m_grid.getDimensions();
    if (i >= dims[0] || j >= dims[1] || k >= dims[2]) {
        throw std::out_of_range("DoseMatrix::at: index out of bounds");
    }
    return m_data[i + j * dims[0] + k * dims[0] * dims[1]];
}

double DoseMatrix::getMax() const {
    if (m_data.empty()) {
        return 0.0;  // Return 0 for empty dose matrix
    }
    return *std::max_element(m_data.begin(), m_data.end());
}

double DoseMatrix::getMean() const {
    if (m_data.empty()) {
        return 0.0;  // Return 0 for empty dose matrix
    }
    return std::accumulate(m_data.begin(), m_data.end(), 0.0) / m_data.size();
}

double DoseMatrix::interpolateAt(double fi, double fj, double fk) const {
    auto dims = m_grid.getDimensions();
    const size_t nx = dims[0], ny = dims[1], nz = dims[2];
    if (nx == 0 || ny == 0 || nz == 0) return 0.0;

    // Clamp to valid range [0, dim-1]
    fi = std::max(0.0, std::min(fi, static_cast<double>(nx - 1)));
    fj = std::max(0.0, std::min(fj, static_cast<double>(ny - 1)));
    fk = std::max(0.0, std::min(fk, static_cast<double>(nz - 1)));

    size_t i0 = static_cast<size_t>(fi);
    size_t j0 = static_cast<size_t>(fj);
    size_t k0 = static_cast<size_t>(fk);
    size_t i1 = std::min(i0 + 1, nx - 1);
    size_t j1 = std::min(j0 + 1, ny - 1);
    size_t k1 = std::min(k0 + 1, nz - 1);

    double di = fi - static_cast<double>(i0);
    double dj = fj - static_cast<double>(j0);
    double dk = fk - static_cast<double>(k0);

    // Flat index: i + j*nx + k*nx*ny
    auto idx = [&](size_t i, size_t j, size_t k) -> size_t {
        return i + j * nx + k * nx * ny;
    };

    double c000 = m_data[idx(i0, j0, k0)];
    double c100 = m_data[idx(i1, j0, k0)];
    double c010 = m_data[idx(i0, j1, k0)];
    double c110 = m_data[idx(i1, j1, k0)];
    double c001 = m_data[idx(i0, j0, k1)];
    double c101 = m_data[idx(i1, j0, k1)];
    double c011 = m_data[idx(i0, j1, k1)];
    double c111 = m_data[idx(i1, j1, k1)];

    // Trilinear interpolation
    double c00 = c000 * (1.0 - di) + c100 * di;
    double c10 = c010 * (1.0 - di) + c110 * di;
    double c01 = c001 * (1.0 - di) + c101 * di;
    double c11 = c011 * (1.0 - di) + c111 * di;

    double c0 = c00 * (1.0 - dj) + c10 * dj;
    double c1 = c01 * (1.0 - dj) + c11 * dj;

    return c0 * (1.0 - dk) + c1 * dk;
}

double* DoseMatrix::data() { return m_data.data(); }
const double* DoseMatrix::data() const { return m_data.data(); }
size_t DoseMatrix::size() const { return m_data.size(); }

} // namespace optirad
