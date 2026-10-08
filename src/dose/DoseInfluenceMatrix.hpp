#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace optirad {

/**
 * Sparse Dose Influence Matrix, stored as a list of CSR column blocks.
 *
 * Construction  – append COO entries (setValue / appendBatch), then close a block with endBlock()
 *                 (the engine closes one block per beam, so only one beam's COO is ever in memory).
 * After finalize() – blocks are read-only; computeDose / accumulateTransposeProduct run over them.
 *
 * Storage is 8 bytes per non-zero (uint32 column + float value). Dose sums are accumulated in double.
 * Row count and column count must each fit in uint32, and each block must hold < 2^32 entries.
 */
class DoseInfluenceMatrix {
public:
    /// One CSR block covering a contiguous or arbitrary set of columns (global column indices).
    struct Block {
        std::vector<uint32_t> rowPtrs;      // numVoxels + 1
        std::vector<uint32_t> colIndices;   // sorted within each row
        std::vector<float>    values;
    };

    DoseInfluenceMatrix() = default;
    DoseInfluenceMatrix(size_t numVoxels, size_t numBixels);

    void setDimensions(size_t numVoxels, size_t numBixels);

    // --- Construction (COO phase) ---
    /// Append a non-zero entry. Duplicate (row,col) are summed when the block is closed.
    void setValue(size_t voxel, size_t bixel, double value);

    /// Bulk-append pre-filtered COO entries (used by parallel ray processing).
    /// Caller is responsible for bounds checking. Not thread-safe; synchronise externally.
    void appendBatch(const std::vector<uint32_t>& rows,
                     const std::vector<uint32_t>& cols,
                     const std::vector<float>& vals);

    /// Convert the pending COO entries into a CSR block and release the COO memory.
    void endBlock();

    /// Reserve COO capacity hint (optional).
    void reserveNonZeros(size_t nnz);

    /// Close any pending block and mark the matrix read-only.
    void finalize();
    bool isFinalized() const { return m_finalized; }

    // --- Access ---
    double getValue(size_t voxel, size_t bixel) const;
    double operator()(size_t voxel, size_t bixel) const;

    size_t getNumVoxels() const;
    size_t getNumBixels() const;
    size_t getNumNonZeros() const;

    // --- Linear algebra (require finalized matrix) ---
    std::vector<double> computeDose(const std::vector<double>& weights) const;
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                    std::vector<double>& grad) const;

    /// Maximum value across all entries (requires finalized matrix).
    double getMaxValue() const;

    // --- Block access (serialization, tests) ---
    const std::vector<Block>& getBlocks() const { return m_blocks; }

    /// Append a pre-built block (deserialization). Marks the matrix as finalized.
    void addBlock(Block block);

    /// Load a single block from CSR arrays. Marks the matrix as finalized.
    void loadCSR(std::vector<size_t> rowPtrs,
                 std::vector<size_t> colIndices,
                 std::vector<double>  values);

    bool isSparse() const { return m_finalized; }

    /// Mark which dose voxels (rows) were actually computed. An empty mask means all rows.
    void setComputedRows(std::vector<uint8_t> mask);
    bool hasAllRows() const { return m_rowMask.empty(); }
    bool isRowComputed(size_t voxel) const { return m_rowMask.empty() || m_rowMask[voxel] != 0; }
    const std::vector<uint8_t>& getRowMask() const { return m_rowMask; }

private:
    size_t m_numVoxels = 0;
    size_t m_numBixels = 0;
    bool   m_finalized = false;

    // COO storage of the block being built
    std::vector<uint32_t> m_cooRows;
    std::vector<uint32_t> m_cooCols;
    std::vector<float>    m_cooVals;

    std::vector<Block> m_blocks;
    std::vector<uint8_t> m_rowMask;  // empty = every row computed
};

} // namespace optirad
