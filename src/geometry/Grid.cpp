#include "Grid.hpp"
#include <stdexcept>
#include <cmath>

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

void Grid::updateMatricesIfNeeded() const {
    if (!m_matricesDirty) return;
    m_matricesDirty = false;
    
    // Row & column direction cosines from DICOM
    m_rowDir = {m_imageOrientation[0], m_imageOrientation[1], m_imageOrientation[2]};
    m_colDir = {m_imageOrientation[3], m_imageOrientation[4], m_imageOrientation[5]};
    
    // Slice direction = cross product
    m_sliceDir = cross(m_rowDir, m_colDir);
    
    // Direction matrix (columns are voxel axes in LPS)
    m_directionMatrix = Mat3(
        m_colDir[0], m_rowDir[0], m_sliceDir[0],
        m_colDir[1], m_rowDir[1], m_sliceDir[1],
        m_colDir[2], m_rowDir[2], m_sliceDir[2]
    );
    
    m_inverseDirectionMatrix = inverse(m_directionMatrix);
}

Vec3 Grid::voxelToPatient(const Vec3& ijk) const {
    updateMatricesIfNeeded();
    
    Vec3 scaled = {
        ijk[0] * m_spacing[0],
        ijk[1] * m_spacing[1],
        ijk[2] * m_spacing[2]
    };
    
    Vec3 rotated = m_directionMatrix * scaled;
    return {
        m_origin[0] + rotated[0],
        m_origin[1] + rotated[1],
        m_origin[2] + rotated[2]
    };
}

Vec3 Grid::patientToVoxel(const Vec3& lps) const {
    updateMatricesIfNeeded();
    
    Vec3 delta = {
        lps[0] - m_origin[0],
        lps[1] - m_origin[1],
        lps[2] - m_origin[2]
    };
    
    Vec3 rotated = m_inverseDirectionMatrix * delta;
    
    // Protect against division by zero
    constexpr double epsilon = 1e-10;
    if (std::abs(m_spacing[0]) < epsilon || 
        std::abs(m_spacing[1]) < epsilon || 
        std::abs(m_spacing[2]) < epsilon) {
        throw std::runtime_error("Grid::patientToVoxel: spacing too small or zero");
    }
    
    return {
        rotated[0] / m_spacing[0],
        rotated[1] / m_spacing[1],
        rotated[2] / m_spacing[2]
    };
}

Vec3 Grid::getRowDirection() const {
    updateMatricesIfNeeded();
    return m_rowDir;
}

Vec3 Grid::getColumnDirection() const {
    updateMatricesIfNeeded();
    return m_colDir;
}

Vec3 Grid::getSliceDirection() const {
    updateMatricesIfNeeded();
    return m_sliceDir;
}

std::vector<double> Grid::getXCoordinates() const {
    std::vector<double> x(m_dimensions[1]);
    for (size_t j = 0; j < m_dimensions[1]; ++j) {
        Vec3 pos = voxelToPatient({0.0, static_cast<double>(j), 0.0});
        x[j] = pos[0];
    }
    return x;
}

std::vector<double> Grid::getYCoordinates() const {
    std::vector<double> y(m_dimensions[0]);
    for (size_t i = 0; i < m_dimensions[0]; ++i) {
        Vec3 pos = voxelToPatient({static_cast<double>(i), 0.0, 0.0});
        y[i] = pos[1];
    }
    return y;
}

std::vector<double> Grid::getZCoordinates() const {
    std::vector<double> z(m_dimensions[2]);
    for (size_t k = 0; k < m_dimensions[2]; ++k) {
        Vec3 pos = voxelToPatient({0.0, 0.0, static_cast<double>(k)});
        z[k] = pos[2];
    }
    return z;
}

} // namespace optirad
