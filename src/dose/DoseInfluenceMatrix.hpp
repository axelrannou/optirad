#pragma once

#include <vector>
#include <cstddef>

namespace optirad {

/**
 * Sparse Dose Influence Matrix (COO build → CSR finalized).
 *
 * Construction phase  – append entries with setValue()  (COO triplets)
 * After finalize()    – CSR for fast computeDose / accumulateTransposeProduct
 *
 * No dense allocation ever – safe for multi-million voxel grids.
 */
class DoseInfluenceMatrix {
public:
    DoseInfluenceMatrix() = default;
    DoseInfluenceMatrix(size_t numVoxels, size_t numBixels);

    void setDimensions(size_t numVoxels, size_t numBixels);

    // --- Construction (COO phase) ---
    /// Append a non-zero entry. Duplicate (row,col) are summed during finalize().
    void setValue(size_t voxel, size_t bixel, double value);

    // --- Access (CSR phase, after finalize()) ---
    double getValue(size_t voxel, size_t bixel) const;
    double operator()(size_t voxel, size_t bixel) const;

    /// Convert COO → CSR. Must be called before computeDose / accumulateTranspose.
    void finalize();
    bool isFinalized() const { return m_finalized; }

    size_t getNumVoxels() const;
    size_t getNumBixels() const;
    size_t getNumNonZeros() const;

    // --- Linear algebra (require finalized CSR) ---
    std::vector<double> computeDose(const std::vector<double>& weights) const;
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                    std::vector<double>& grad) const;

    // --- CSR direct access (for serialization) ---
    const std::vector<double>&  getValues()     const { return m_values; }
    const std::vector<size_t>&  getColIndices() const { return m_colIndices; }
    const std::vector<size_t>&  getRowPtrs()    const { return m_rowPtrs; }

    /// Load pre-built CSR arrays directly (deserialization). Marks matrix as finalized.
    void loadCSR(std::vector<size_t> rowPtrs,
                 std::vector<size_t> colIndices,
                 std::vector<double>  values);

    bool isSparse() const { return m_finalized; }

    /// Reserve COO capacity hint (optional, avoids reallocation when nnz is known).
    void reserveNonZeros(size_t nnz);

    /// Bulk-append pre-filtered COO entries (used by parallel ray processing).
    /// Caller is responsible for bounds checking.  Not thread-safe on its own;
    /// use external synchronisation (e.g., omp critical) when called from threads.
    void appendBatch(const std::vector<size_t>& rows,
                     const std::vector<size_t>& cols,
                     const std::vector<double>& vals);

private:
    size_t m_numVoxels = 0;
    size_t m_numBixels = 0;
    bool   m_finalized = false;

    // COO storage (during construction)
    std::vector<size_t> m_cooRows;
    std::vector<size_t> m_cooCols;
    std::vector<double> m_cooVals;

    // CSR storage (after finalize)
    std::vector<double>  m_values;
    std::vector<size_t>  m_colIndices;
    std::vector<size_t>  m_rowPtrs;
};

} // namespace optirad
