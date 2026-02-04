#include "DoseInfluenceMatrix.hpp"
#include <algorithm>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

DoseInfluenceMatrix::DoseInfluenceMatrix(size_t numVoxels, size_t numBixels) {
    setDimensions(numVoxels, numBixels);
}

void DoseInfluenceMatrix::setDimensions(size_t numVoxels, size_t numBixels) {
    m_numVoxels = numVoxels;
    m_numBixels = numBixels;
    m_denseData.resize(numVoxels * numBixels, 0.0);
    m_useSparse = false;
}

double& DoseInfluenceMatrix::operator()(size_t voxel, size_t bixel) {
    return m_denseData[voxel * m_numBixels + bixel];
}

double DoseInfluenceMatrix::operator()(size_t voxel, size_t bixel) const {
    return getValue(voxel, bixel);
}

void DoseInfluenceMatrix::setValue(size_t voxel, size_t bixel, double value) {
    if (m_useSparse) {
        throw std::runtime_error("Cannot setValue after finalize()");
    }
    m_denseData[voxel * m_numBixels + bixel] = value;
}

double DoseInfluenceMatrix::getValue(size_t voxel, size_t bixel) const {
    if (!m_useSparse) {
        return m_denseData[voxel * m_numBixels + bixel];
    }
    // Binary search in sparse row
    size_t rowStart = m_rowPtrs[voxel];
    size_t rowEnd = m_rowPtrs[voxel + 1];
    auto it = std::lower_bound(m_colIndices.begin() + rowStart,
                               m_colIndices.begin() + rowEnd, bixel);
    if (it != m_colIndices.begin() + rowEnd && *it == bixel) {
        return m_values[it - m_colIndices.begin()];
    }
    return 0.0;
}

void DoseInfluenceMatrix::finalize() {
    // Count non-zeros
    size_t nnz = 0;
    for (double v : m_denseData) {
        if (v != 0.0) ++nnz;
    }
    
    double sparsity = 1.0 - static_cast<double>(nnz) / m_denseData.size();
    
    // Only convert to sparse if beneficial
    if (sparsity > (1.0 - SPARSE_THRESHOLD) && m_denseData.size() > 10000) {
        m_values.reserve(nnz);
        m_colIndices.reserve(nnz);
        m_rowPtrs.resize(m_numVoxels + 1);
        
        size_t idx = 0;
        for (size_t v = 0; v < m_numVoxels; ++v) {
            m_rowPtrs[v] = m_values.size();
            for (size_t b = 0; b < m_numBixels; ++b) {
                double val = m_denseData[v * m_numBixels + b];
                if (val != 0.0) {
                    m_values.push_back(val);
                    m_colIndices.push_back(b);
                }
            }
        }
        m_rowPtrs[m_numVoxels] = m_values.size();
        
        m_denseData.clear();
        m_denseData.shrink_to_fit();
        m_useSparse = true;
    }
}

size_t DoseInfluenceMatrix::getNumVoxels() const { return m_numVoxels; }
size_t DoseInfluenceMatrix::getNumBixels() const { return m_numBixels; }
size_t DoseInfluenceMatrix::getNumNonZeros() const {
    return m_useSparse ? m_values.size() : 
           std::count_if(m_denseData.begin(), m_denseData.end(), 
                        [](double v) { return v != 0.0; });
}

std::vector<double> DoseInfluenceMatrix::computeDose(const std::vector<double>& weights) const {
    std::vector<double> dose(m_numVoxels, 0.0);
    
    if (m_useSparse) {
        #pragma omp parallel for
        for (size_t v = 0; v < m_numVoxels; ++v) {
            double sum = 0.0;
            for (size_t k = m_rowPtrs[v]; k < m_rowPtrs[v + 1]; ++k) {
                sum += m_values[k] * weights[m_colIndices[k]];
            }
            dose[v] = sum;
        }
    } else {
        #pragma omp parallel for
        for (size_t v = 0; v < m_numVoxels; ++v) {
            double sum = 0.0;
            for (size_t b = 0; b < m_numBixels; ++b) {
                sum += m_denseData[v * m_numBixels + b] * weights[b];
            }
            dose[v] = sum;
        }
    }
    return dose;
}

void DoseInfluenceMatrix::accumulateTransposeProduct(
    const std::vector<double>& voxelGrad,
    std::vector<double>& grad) const 
{
    if (m_useSparse) {
        // For sparse: iterate by row (voxel), accumulate to columns (beamlets)
        // Need thread-safe accumulation
        #pragma omp parallel for
        for (size_t v = 0; v < m_numVoxels; ++v) {
            double gv = voxelGrad[v];
            if (gv == 0.0) continue;
            for (size_t k = m_rowPtrs[v]; k < m_rowPtrs[v + 1]; ++k) {
                #pragma omp atomic
                grad[m_colIndices[k]] += gv * m_values[k];
            }
        }
    } else {
        #pragma omp parallel for
        for (size_t b = 0; b < m_numBixels; ++b) {
            double sum = 0.0;
            for (size_t v = 0; v < m_numVoxels; ++v) {
                sum += voxelGrad[v] * m_denseData[v * m_numBixels + b];
            }
            grad[b] += sum;
        }
    }
}

} // namespace optirad
