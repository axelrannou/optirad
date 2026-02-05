#pragma once

#include "MathUtils.hpp"
#include <array>
#include <cstddef>
#include <string>

namespace optirad {

class Grid {
public:
    void setDimensions(size_t nx, size_t ny, size_t nz);
    void setSpacing(double dx, double dy, double dz);
    void setOrigin(const Vec3& origin);
    
    // DICOM geometric info
    void setPatientPosition(const std::string& pos) { m_patientPosition = pos; }
    void setImageOrientation(const std::array<double, 6>& orient) { m_imageOrientation = orient; }
    void setSliceThickness(double thickness) { m_sliceThickness = thickness; }

    std::array<size_t, 3> getDimensions() const;
    Vec3 getSpacing() const;
    Vec3 getOrigin() const;
    size_t getNumVoxels() const;
    
    const std::string& getPatientPosition() const { return m_patientPosition; }
    const std::array<double, 6>& getImageOrientation() const { return m_imageOrientation; }
    double getSliceThickness() const { return m_sliceThickness; }

private:
    std::array<size_t, 3> m_dimensions = {0, 0, 0};
    Vec3 m_spacing = {1.0, 1.0, 1.0};
    Vec3 m_origin = {0.0, 0.0, 0.0};
    
    // DICOM patient coordinate system
    std::string m_patientPosition = "HFS";  // Head First Supine (default)
    std::array<double, 6> m_imageOrientation = {1, 0, 0, 0, 1, 0};  // Row, Column direction cosines
    double m_sliceThickness = 0.0;  // May differ from spacing[2]
};

} // namespace optirad
