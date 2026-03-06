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

Grid Grid::createDoseGrid(const Grid& ctGrid, const Vec3& doseResolution) {
    auto ctDims = ctGrid.getDimensions();
    auto ctSpacing = ctGrid.getSpacing();
    auto ctOrigin = ctGrid.getOrigin();

    // Match matRad behavior:
    // - anchor dose grid at first CT voxel
    // - sample with step = dose resolution up to CT max extent
    // This corresponds to x = ct.x(1):dx:ct.x(end), same for y/z.
    // Therefore: N = floor((ctMax - ctMin)/doseRes) + 1.
    double extentX = static_cast<double>(ctDims[0] - 1) * ctSpacing[0];
    double extentY = static_cast<double>(ctDims[1] - 1) * ctSpacing[1];
    double extentZ = static_cast<double>(ctDims[2] - 1) * ctSpacing[2];

    size_t doseNx = std::max(size_t(1), static_cast<size_t>(std::floor(extentX / doseResolution[0])) + 1);
    size_t doseNy = std::max(size_t(1), static_cast<size_t>(std::floor(extentY / doseResolution[1])) + 1);
    size_t doseNz = std::max(size_t(1), static_cast<size_t>(std::floor(extentZ / doseResolution[2])) + 1);

    Grid doseGrid;
    doseGrid.setDimensions(doseNx, doseNy, doseNz);
    doseGrid.setSpacing(doseResolution[0], doseResolution[1], doseResolution[2]);
    doseGrid.setPatientPosition(ctGrid.getPatientPosition());
    doseGrid.setImageOrientation(ctGrid.getImageOrientation());
    doseGrid.setSliceThickness(doseResolution[2]);

    // Anchor to CT origin (first voxel world coordinate)
    doseGrid.setOrigin(ctOrigin);

    return doseGrid;
}

std::vector<size_t> Grid::mapVoxelIndices(const Grid& fromGrid, const Grid& toGrid,
                                           const std::vector<size_t>& fromIndices) {
    auto fromDims = fromGrid.getDimensions();
    auto toDims = toGrid.getDimensions();
    size_t fromNx = fromDims[0], fromNy = fromDims[1];
    size_t toNx = toDims[0], toNy = toDims[1], toNz = toDims[2];

    std::vector<size_t> mapped;
    mapped.reserve(fromIndices.size());

    // Use a set to avoid duplicates
    std::vector<bool> seen(toGrid.getNumVoxels(), false);

    for (size_t idx : fromIndices) {
        // Convert flat index to ijk in source grid
        size_t iz = idx / (fromNx * fromNy);
        size_t rem = idx % (fromNx * fromNy);
        size_t iy = rem / fromNx;
        size_t ix = rem % fromNx;

        // Convert to patient coordinates
        Vec3 lps = fromGrid.voxelToPatient({
            static_cast<double>(ix),
            static_cast<double>(iy),
            static_cast<double>(iz)
        });

        // Convert to target grid voxel coordinates
        Vec3 ijk = toGrid.patientToVoxel(lps);
        int tix = static_cast<int>(std::round(ijk[0]));
        int tiy = static_cast<int>(std::round(ijk[1]));
        int tiz = static_cast<int>(std::round(ijk[2]));

        // Bounds check
        if (tix < 0 || tix >= static_cast<int>(toNx) ||
            tiy < 0 || tiy >= static_cast<int>(toNy) ||
            tiz < 0 || tiz >= static_cast<int>(toNz)) {
            continue;
        }

        size_t toIdx = static_cast<size_t>(tix) +
                       static_cast<size_t>(tiy) * toNx +
                       static_cast<size_t>(tiz) * toNx * toNy;

        if (!seen[toIdx]) {
            seen[toIdx] = true;
            mapped.push_back(toIdx);
        }
    }

    return mapped;
}

std::vector<size_t> Grid::resampleMaskNearestToGrid(const Grid& sourceGrid, const Grid& targetGrid,
                                                    const std::vector<size_t>& sourceMaskIndices) {
    auto sourceDims = sourceGrid.getDimensions();
    auto targetDims = targetGrid.getDimensions();

    const size_t sourceNx = sourceDims[0];
    const size_t sourceNy = sourceDims[1];
    const size_t sourceNz = sourceDims[2];

    const size_t targetNx = targetDims[0];
    const size_t targetNy = targetDims[1];
    const size_t targetNz = targetDims[2];

    std::vector<unsigned char> sourceMask(sourceGrid.getNumVoxels(), 0);
    for (size_t idx : sourceMaskIndices) {
        if (idx < sourceMask.size()) {
            sourceMask[idx] = 1;
        }
    }

    std::vector<size_t> mapped;
    mapped.reserve(targetGrid.getNumVoxels() / 8);

    for (size_t tz = 0; tz < targetNz; ++tz) {
        for (size_t ty = 0; ty < targetNy; ++ty) {
            for (size_t tx = 0; tx < targetNx; ++tx) {
                const size_t targetIdx = tx + ty * targetNx + tz * targetNx * targetNy;

                Vec3 lps = targetGrid.voxelToPatient({
                    static_cast<double>(tx),
                    static_cast<double>(ty),
                    static_cast<double>(tz)
                });

                Vec3 sourceIJK = sourceGrid.patientToVoxel(lps);
                int sx = static_cast<int>(std::round(sourceIJK[0]));
                int sy = static_cast<int>(std::round(sourceIJK[1]));
                int sz = static_cast<int>(std::round(sourceIJK[2]));

                if (sx < 0 || sx >= static_cast<int>(sourceNx) ||
                    sy < 0 || sy >= static_cast<int>(sourceNy) ||
                    sz < 0 || sz >= static_cast<int>(sourceNz)) {
                    continue;
                }

                const size_t sourceIdx = static_cast<size_t>(sx) +
                                         static_cast<size_t>(sy) * sourceNx +
                                         static_cast<size_t>(sz) * sourceNx * sourceNy;

                if (sourceMask[sourceIdx]) {
                    mapped.push_back(targetIdx);
                }
            }
        }
    }

    return mapped;
}

} // namespace optirad
