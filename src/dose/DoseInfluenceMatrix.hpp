#pragma once

#include <vector>
#include <cstddef>

namespace optirad {

/**
 * Sparse Dose Influence Matrix using CSR (Compressed Sparse Row) format
 * Stores D_ij = dose to voxel i from unit fluence of beamlet j
 */
class DoseInfluenceMatrix {
public:
    void setDimensions(size_t numVoxels, size_t numBixels);
    
    // For building the matrix (dense interface - converts to sparse internally)
    void setValue(size_t voxel, size_t bixel, double value);
    double getValue(size_t voxel, size_t bixel) const;
    
    // Finalize sparse structure after all setValue calls
    void finalize();
    
    size_t getNumVoxels() const;
    size_t getNumBixels() const;
    size_t getNumNonZeros() const;
    
    // Compute dose = Dij * weights (optimized for sparse)
    std::vector<double> computeDose(const std::vector<double>& weights) const;
    
    // Compute gradient contribution: grad += Dij^T * voxelGrad
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                    std::vector<double>& grad) const;
    
    // Check if using sparse format
    bool isSparse() const { return m_useSparse; }

private:
    size_t m_numVoxels = 0;
    size_t m_numBixels = 0;
    
    // Dense storage (used during construction or for small matrices)
    std::vector<double> m_denseData;
    
    // CSR sparse storage
    bool m_useSparse = false;
    std::vector<double> m_values;      // Non-zero values
    std::vector<size_t> m_colIndices;  // Column index for each value
    std::vector<size_t> m_rowPtrs;     // Start of each row in values/colIndices
    
    static constexpr double SPARSE_THRESHOLD = 0.1;  // Use sparse if <10% non-zeros
};

} // namespace optirad
