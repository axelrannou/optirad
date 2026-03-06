#include "DoseInfluenceMatrix.hpp"
#include <algorithm>
#include <stdexcept>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

// ────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────

DoseInfluenceMatrix::DoseInfluenceMatrix(size_t numVoxels, size_t numBixels)
    : m_numVoxels(numVoxels), m_numBixels(numBixels) {}

void DoseInfluenceMatrix::setDimensions(size_t numVoxels, size_t numBixels) {
    m_numVoxels = numVoxels;
    m_numBixels = numBixels;
    m_finalized = false;
    m_cooRows.clear();
    m_cooCols.clear();
    m_cooVals.clear();
    m_values.clear();
    m_colIndices.clear();
    m_rowPtrs.clear();
}

void DoseInfluenceMatrix::reserveNonZeros(size_t nnz) {
    m_cooRows.reserve(nnz);
    m_cooCols.reserve(nnz);
    m_cooVals.reserve(nnz);
}

void DoseInfluenceMatrix::appendBatch(const std::vector<size_t>& rows,
                                       const std::vector<size_t>& cols,
                                       const std::vector<double>& vals) {
    if (m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix: cannot appendBatch after finalize()");
    m_cooRows.insert(m_cooRows.end(), rows.begin(), rows.end());
    m_cooCols.insert(m_cooCols.end(), cols.begin(), cols.end());
    m_cooVals.insert(m_cooVals.end(), vals.begin(), vals.end());
}

// ────────────────────────────────────────────────────────────────
// COO accumulation
// ────────────────────────────────────────────────────────────────

void DoseInfluenceMatrix::setValue(size_t voxel, size_t bixel, double value) {
    if (m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix: cannot setValue after finalize()");
    if (voxel >= m_numVoxels || bixel >= m_numBixels)
        throw std::out_of_range("DoseInfluenceMatrix::setValue: index out of bounds");
    m_cooRows.push_back(voxel);
    m_cooCols.push_back(bixel);
    m_cooVals.push_back(value);
}

// ────────────────────────────────────────────────────────────────
// Finalize: COO → CSR
// ────────────────────────────────────────────────────────────────

void DoseInfluenceMatrix::finalize() {
    if (m_finalized) return;

    size_t nnz = m_cooRows.size();

    // Sort COO by (row, col) using an index permutation
    std::vector<size_t> perm(nnz);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b) {
        if (m_cooRows[a] != m_cooRows[b]) return m_cooRows[a] < m_cooRows[b];
        return m_cooCols[a] < m_cooCols[b];
    });

    // Build CSR, merging duplicate (row,col) entries by summing values
    m_rowPtrs.assign(m_numVoxels + 1, 0);
    m_values.clear();
    m_colIndices.clear();
    m_values.reserve(nnz);
    m_colIndices.reserve(nnz);

    for (size_t pi = 0; pi < nnz; ++pi) {
        size_t i = perm[pi];
        size_t row = m_cooRows[i];
        size_t col = m_cooCols[i];
        double val = m_cooVals[i];

        // Merge with previous entry if same (row, col)
        if (!m_values.empty() && !m_colIndices.empty() &&
            m_colIndices.back() == col &&
            m_rowPtrs[row + 1] == 0 &&  // still on same row check below
            pi > 0 && m_cooRows[perm[pi - 1]] == row)
        {
            // Simpler merge check: if the last pushed entry has same col and is in same row
            m_values.back() += val;
        } else {
            m_values.push_back(val);
            m_colIndices.push_back(col);
        }
    }

    // Actually, the merge logic above is fragile. Let me redo it cleanly.
    // Rebuild from scratch with a cleaner approach:
    m_values.clear();
    m_colIndices.clear();
    m_rowPtrs.assign(m_numVoxels + 1, 0);

    if (nnz > 0) {
        size_t prevRow = m_cooRows[perm[0]];
        size_t prevCol = m_cooCols[perm[0]];
        double accum  = m_cooVals[perm[0]];

        for (size_t pi = 1; pi < nnz; ++pi) {
            size_t i = perm[pi];
            size_t row = m_cooRows[i];
            size_t col = m_cooCols[i];

            if (row == prevRow && col == prevCol) {
                accum += m_cooVals[i]; // merge duplicate
            } else {
                // Flush previous entry
                m_values.push_back(accum);
                m_colIndices.push_back(prevCol);
                // Mark row boundaries
                for (size_t r = prevRow + 1; r <= row; ++r) {
                    // rows [prevRow+1 .. row] start at current position
                }
                prevRow = row;
                prevCol = col;
                accum = m_cooVals[i];
            }
        }
        // Flush last entry
        m_values.push_back(accum);
        m_colIndices.push_back(prevCol);
    }

    // Build row pointers by counting entries per row
    // (redo with a simple counting approach)
    m_rowPtrs.assign(m_numVoxels + 1, 0);
    {
        // Count entries per row from the deduplicated sorted COO
        // We have m_values and m_colIndices, but we need to reconstruct row info.
        // Let me redo completely with a cleaner 2-pass approach.
    }

    // --- CLEAN 2-pass COO → CSR ---
    m_values.clear();
    m_colIndices.clear();
    m_rowPtrs.assign(m_numVoxels + 1, 0);

    // Pass 1: count entries per row (after dedup)
    if (nnz > 0) {
        size_t prevRow = m_cooRows[perm[0]];
        size_t prevCol = m_cooCols[perm[0]];
        for (size_t pi = 1; pi < nnz; ++pi) {
            size_t i = perm[pi];
            if (m_cooRows[i] != prevRow || m_cooCols[i] != prevCol) {
                m_rowPtrs[prevRow + 1]++; // count for row prevRow
                prevRow = m_cooRows[i];
                prevCol = m_cooCols[i];
            }
        }
        m_rowPtrs[prevRow + 1]++; // last entry
    }

    // Cumulative sum
    for (size_t r = 0; r < m_numVoxels; ++r) {
        m_rowPtrs[r + 1] += m_rowPtrs[r];
    }

    size_t totalNnz = m_rowPtrs[m_numVoxels];
    m_values.resize(totalNnz);
    m_colIndices.resize(totalNnz);

    // Pass 2: fill values and column indices
    if (nnz > 0) {
        std::vector<size_t> rowCursor(m_rowPtrs.begin(), m_rowPtrs.begin() + m_numVoxels);

        size_t prevRow = m_cooRows[perm[0]];
        size_t prevCol = m_cooCols[perm[0]];
        double accum = m_cooVals[perm[0]];

        for (size_t pi = 1; pi < nnz; ++pi) {
            size_t i = perm[pi];
            size_t row = m_cooRows[i];
            size_t col = m_cooCols[i];

            if (row == prevRow && col == prevCol) {
                accum += m_cooVals[i];
            } else {
                size_t pos = rowCursor[prevRow]++;
                m_values[pos] = accum;
                m_colIndices[pos] = prevCol;
                prevRow = row;
                prevCol = col;
                accum = m_cooVals[i];
            }
        }
        size_t pos = rowCursor[prevRow]++;
        m_values[pos] = accum;
        m_colIndices[pos] = prevCol;
    }

    // Free COO memory
    m_cooRows.clear();   m_cooRows.shrink_to_fit();
    m_cooCols.clear();   m_cooCols.shrink_to_fit();
    m_cooVals.clear();   m_cooVals.shrink_to_fit();

    m_finalized = true;
}

// ────────────────────────────────────────────────────────────────
// Direct CSR loading (deserialization)
// ────────────────────────────────────────────────────────────────

void DoseInfluenceMatrix::loadCSR(std::vector<size_t> rowPtrs,
                                   std::vector<size_t> colIndices,
                                   std::vector<double>  values) {
    m_rowPtrs    = std::move(rowPtrs);
    m_colIndices = std::move(colIndices);
    m_values     = std::move(values);
    m_finalized  = true;
    // Free any leftover COO
    m_cooRows.clear(); m_cooRows.shrink_to_fit();
    m_cooCols.clear(); m_cooCols.shrink_to_fit();
    m_cooVals.clear(); m_cooVals.shrink_to_fit();
}

// ────────────────────────────────────────────────────────────────
// Read-only access
// ────────────────────────────────────────────────────────────────

double DoseInfluenceMatrix::operator()(size_t voxel, size_t bixel) const {
    return getValue(voxel, bixel);
}

double DoseInfluenceMatrix::getValue(size_t voxel, size_t bixel) const {
    if (voxel >= m_numVoxels || bixel >= m_numBixels)
        throw std::out_of_range("DoseInfluenceMatrix::getValue: index out of bounds");

    if (!m_finalized) {
        // Linear scan through COO (slow – only for debugging/tests)
        double sum = 0.0;
        for (size_t k = 0; k < m_cooRows.size(); ++k) {
            if (m_cooRows[k] == voxel && m_cooCols[k] == bixel)
                sum += m_cooVals[k];
        }
        return sum;
    }

    // Binary search in CSR row
    size_t rowStart = m_rowPtrs[voxel];
    size_t rowEnd   = m_rowPtrs[voxel + 1];
    auto beg = m_colIndices.begin() + static_cast<ptrdiff_t>(rowStart);
    auto end = m_colIndices.begin() + static_cast<ptrdiff_t>(rowEnd);
    auto it  = std::lower_bound(beg, end, bixel);
    if (it != end && *it == bixel)
        return m_values[static_cast<size_t>(it - m_colIndices.begin())];
    return 0.0;
}

// ────────────────────────────────────────────────────────────────
// Dimensions / stats
// ────────────────────────────────────────────────────────────────

size_t DoseInfluenceMatrix::getNumVoxels() const { return m_numVoxels; }
size_t DoseInfluenceMatrix::getNumBixels() const { return m_numBixels; }

size_t DoseInfluenceMatrix::getNumNonZeros() const {
    return m_finalized ? m_values.size() : m_cooRows.size();
}

double DoseInfluenceMatrix::getMaxValue() const {
    if (!m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix::getMaxValue requires finalize()");
    if (m_values.empty()) return 0.0;
    return *std::max_element(m_values.begin(), m_values.end());
}

// ────────────────────────────────────────────────────────────────
// Linear algebra (all CSR-based, require finalized)
// ────────────────────────────────────────────────────────────────

std::vector<double> DoseInfluenceMatrix::computeDose(const std::vector<double>& weights) const {
    if (!m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix::computeDose requires finalize()");

    std::vector<double> dose(m_numVoxels, 0.0);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t v = 0; v < m_numVoxels; ++v) {
        double sum = 0.0;
        for (size_t k = m_rowPtrs[v]; k < m_rowPtrs[v + 1]; ++k)
            sum += m_values[k] * weights[m_colIndices[k]];
        dose[v] = sum;
    }
    return dose;
}

void DoseInfluenceMatrix::accumulateTransposeProduct(
    const std::vector<double>& voxelGrad,
    std::vector<double>& grad) const
{
    if (!m_finalized)
        throw std::runtime_error("accumulateTransposeProduct requires finalize()");

    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t v = 0; v < m_numVoxels; ++v) {
        double gv = voxelGrad[v];
        if (gv == 0.0) continue;
        for (size_t k = m_rowPtrs[v]; k < m_rowPtrs[v + 1]; ++k) {
            #pragma omp atomic
            grad[m_colIndices[k]] += gv * m_values[k];
        }
    }
}

} // namespace optirad
