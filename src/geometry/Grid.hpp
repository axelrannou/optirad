#pragma once

#include "MathUtils.hpp"  // Local, not utils/MathUtils.hpp
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace optirad {

class Grid {
public:
    void setDimensions(size_t nx, size_t ny, size_t nz);
    void setSpacing(double dx, double dy, double dz);
    void setOrigin(const Vec3& origin);
    
    // DICOM geometric info
    void setPatientPosition(const std::string& pos) { m_patientPosition = pos; }
    void setImageOrientation(const std::array<double, 6>& orient) { 
        m_imageOrientation = orient; 
        m_matricesDirty = true;
    }
    void setSliceThickness(double thickness) { m_sliceThickness = thickness; }

    std::array<size_t, 3> getDimensions() const;
    Vec3 getSpacing() const;
    Vec3 getOrigin() const;
    size_t getNumVoxels() const;
    
    const std::string& getPatientPosition() const { return m_patientPosition; }
    const std::array<double, 6>& getImageOrientation() const { return m_imageOrientation; }
    double getSliceThickness() const { return m_sliceThickness; }
    
    // Coordinate transforms
    Vec3 voxelToPatient(const Vec3& ijk) const;
    Vec3 patientToVoxel(const Vec3& lps) const;
    
    // Direction accessors
    Vec3 getRowDirection() const;
    Vec3 getColumnDirection() const;
    Vec3 getSliceDirection() const;
    
    // Coordinate arrays
    std::vector<double> getXCoordinates() const;
    std::vector<double> getYCoordinates() const;
    std::vector<double> getZCoordinates() const;

    /// Create a dose calculation grid from a CT grid with different resolution.
    /// The dose grid covers the same physical extent but with coarser spacing.
    static Grid createDoseGrid(const Grid& ctGrid, const Vec3& doseResolution);

    /// Map flat voxel indices from one grid to another via nearest-neighbor.
    /// Returns dose-grid flat indices corresponding to ctVoxelIndices.
    static std::vector<size_t> mapVoxelIndices(const Grid& fromGrid, const Grid& toGrid,
                                                const std::vector<size_t>& fromIndices);

private:
    void updateMatricesIfNeeded() const;
    
    std::array<size_t, 3> m_dimensions = {0, 0, 0};
    Vec3 m_spacing = {1.0, 1.0, 1.0};
    Vec3 m_origin = {0.0, 0.0, 0.0};
    
    // DICOM patient coordinate system
    std::string m_patientPosition = "HFS";
    std::array<double, 6> m_imageOrientation = {1, 0, 0, 0, 1, 0};
    double m_sliceThickness = 0.0;
    
    // Direction cosines (cached)
    mutable Vec3 m_rowDir{1, 0, 0};
    mutable Vec3 m_colDir{0, 1, 0};
    mutable Vec3 m_sliceDir{0, 0, 1};
    
    // Cached matrices
    mutable Mat3 m_directionMatrix;
    mutable Mat3 m_inverseDirectionMatrix;
    mutable bool m_matricesDirty = true;
};

} // namespace optirad
