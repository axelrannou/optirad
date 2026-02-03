#include "DoseInfluenceMatrix.hpp"

namespace optirad {

void DoseInfluenceMatrix::setDimensions(size_t numVoxels, size_t numBixels) {
    m_numVoxels = numVoxels;
    m_numBixels = numBixels;
    m_data.resize(numVoxels * numBixels, 0.0);
}

void DoseInfluenceMatrix::setValue(size_t voxel, size_t bixel, double value) {
    m_data[voxel * m_numBixels + bixel] = value;
}

double DoseInfluenceMatrix::getValue(size_t voxel, size_t bixel) const {
    return m_data[voxel * m_numBixels + bixel];
}

size_t DoseInfluenceMatrix::getNumVoxels() const { return m_numVoxels; }
size_t DoseInfluenceMatrix::getNumBixels() const { return m_numBixels; }

std::vector<double> DoseInfluenceMatrix::computeDose(const std::vector<double>& weights) const {
    std::vector<double> dose(m_numVoxels, 0.0);
    for (size_t v = 0; v < m_numVoxels; ++v) {
        for (size_t b = 0; b < m_numBixels; ++b) {
            dose[v] += getValue(v, b) * weights[b];
        }
    }
    return dose;
}

} // namespace optirad
