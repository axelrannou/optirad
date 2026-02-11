#include "DoseMatrix.hpp"
#include <algorithm>
#include <numeric>

namespace optirad {

void DoseMatrix::setGrid(const Grid& grid) { m_grid = grid; }

void DoseMatrix::allocate() {
    m_data.resize(m_grid.getNumVoxels(), 0.0);
}

double& DoseMatrix::at(size_t i, size_t j, size_t k) {
    auto dims = m_grid.getDimensions();
    return m_data[i + dims[0] * (j + dims[1] * k)];
}

const double& DoseMatrix::at(size_t i, size_t j, size_t k) const {
    auto dims = m_grid.getDimensions();
    return m_data[i + dims[0] * (j + dims[1] * k)];
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

double* DoseMatrix::data() { return m_data.data(); }
const double* DoseMatrix::data() const { return m_data.data(); }
size_t DoseMatrix::size() const { return m_data.size(); }

} // namespace optirad
