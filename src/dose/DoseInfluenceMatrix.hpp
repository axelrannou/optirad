#pragma once

#include <vector>
#include <cstddef>

namespace optirad {

// Sparse dose influence matrix: dose = Dij * weights
// Dij[voxel][bixel] = dose contribution to voxel from unit weight bixel
class DoseInfluenceMatrix {
public:
    void setDimensions(size_t numVoxels, size_t numBixels);
    
    void setValue(size_t voxel, size_t bixel, double value);
    double getValue(size_t voxel, size_t bixel) const;

    size_t getNumVoxels() const;
    size_t getNumBixels() const;

    // Compute dose = Dij * weights
    std::vector<double> computeDose(const std::vector<double>& weights) const;

private:
    size_t m_numVoxels = 0;
    size_t m_numBixels = 0;
    std::vector<double> m_data; // Dense for now, sparse later
};

} // namespace optirad
