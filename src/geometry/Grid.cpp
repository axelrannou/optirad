#include "Grid.hpp"

namespace optirad {

void Grid::setDimensions(size_t nx, size_t ny, size_t nz) {
    m_dimensions = {nx, ny, nz};
}

void Grid::setSpacing(double dx, double dy, double dz) {
    m_spacing = {dx, dy, dz};
}

void Grid::setOrigin(const Vec3& origin) { m_origin = origin; }

std::array<size_t, 3> Grid::getDimensions() const { return m_dimensions; }
Vec3 Grid::getSpacing() const { return m_spacing; }
Vec3 Grid::getOrigin() const { return m_origin; }

size_t Grid::getNumVoxels() const {
    return m_dimensions[0] * m_dimensions[1] * m_dimensions[2];
}

} // namespace optirad
