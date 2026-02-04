#pragma once

#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace optirad {

/**
 * A single contour on a slice (list of 3D points)
 */
struct Contour {
    std::vector<std::array<double, 3>> points;  // x, y, z coordinates
    double zPosition = 0.0;  // Z position of this contour slice
};

/**
 * Structure/ROI definition with contours and voxel mask
 */
class Structure {
public:
    Structure() = default;
    explicit Structure(const std::string& name) : m_name(name) {}
    
    // Basic properties
    void setName(const std::string& name) { m_name = name; }
    const std::string& getName() const { return m_name; }
    
    void setType(const std::string& type) { m_type = type; }
    const std::string& getType() const { return m_type; }
    
    void setROINumber(int num) { m_roiNumber = num; }
    int getROINumber() const { return m_roiNumber; }
    
    void setColor(uint8_t r, uint8_t g, uint8_t b) { m_color = {r, g, b}; }
    std::array<uint8_t, 3> getColor() const { return m_color; }
    
    // Contour data
    void addContour(const Contour& contour) { m_contours.push_back(contour); }
    const std::vector<Contour>& getContours() const { return m_contours; }
    size_t getContourCount() const { return m_contours.size(); }
    
    // Voxel indices (computed from contours, like matRad's cst{i,4})
    void setVoxelIndices(const std::vector<size_t>& indices) { m_voxelIndices = indices; }
    const std::vector<size_t>& getVoxelIndices() const { return m_voxelIndices; }
    size_t getVoxelCount() const { return m_voxelIndices.size(); }
    
    // Optimization parameters (like matRad's cst{i,5})
    void setPriority(int priority) { m_priority = priority; }
    int getPriority() const { return m_priority; }
    
    void setAlphaX(double alpha) { m_alphaX = alpha; }
    double getAlphaX() const { return m_alphaX; }
    
    void setBetaX(double beta) { m_betaX = beta; }
    double getBetaX() const { return m_betaX; }
    
    void setVisible(bool visible) { m_visible = visible; }
    bool isVisible() const { return m_visible; }
    
private:
    std::string m_name;
    std::string m_type = "UNKNOWN";  // TARGET, OAR, EXTERNAL, UNKNOWN
    int m_roiNumber = 0;
    std::array<uint8_t, 3> m_color = {255, 0, 0};
    
    std::vector<Contour> m_contours;
    std::vector<size_t> m_voxelIndices;
    
    // Optimization parameters
    int m_priority = 5;
    double m_alphaX = 0.1;
    double m_betaX = 0.05;
    bool m_visible = true;
};

} // namespace optirad
