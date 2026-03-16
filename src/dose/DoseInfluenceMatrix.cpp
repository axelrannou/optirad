#include "DoseInfluenceMatrix.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

namespace {

using ColVal = std::pair<size_t, double>;

} // namespace

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

    const size_t nnz = m_cooRows.size();
    const size_t cooBytes = nnz * (2 * sizeof(size_t) + sizeof(double));
    const size_t bucketBytes = nnz * sizeof(ColVal);
    const size_t estimatedGB = (cooBytes + bucketBytes) / (1024ull * 1024ull * 1024ull);

    Logger::info("DoseInfluenceMatrix::finalize: Converting " + std::to_string(nnz) +
        " entries to CSR using row bucketing (estimated peak memory: ~" +
        std::to_string(estimatedGB) + " GB before CSR allocation)");

    if (nnz == 0) {
        Logger::warn("DoseInfluenceMatrix::finalize: No entries to finalize");
        m_rowPtrs.assign(m_numVoxels + 1, 0);
        m_finalized = true;
        return;
    }

    Logger::info("DoseInfluenceMatrix::finalize: Counting entries per row...");
    m_rowPtrs.assign(m_numVoxels + 1, 0);
    for (size_t i = 0; i < nnz; ++i) {
        const size_t row = m_cooRows[i];
        const size_t col = m_cooCols[i];
        if (row >= m_numVoxels || col >= m_numBixels) {
            throw std::runtime_error(
                "DoseInfluenceMatrix::finalize: index out of bounds at entry " +
                std::to_string(i) + " (row=" + std::to_string(row) +
                ", col=" + std::to_string(col) + ")");
        }
        ++m_rowPtrs[row + 1];
    }

    for (size_t row = 0; row < m_numVoxels; ++row) {
        m_rowPtrs[row + 1] += m_rowPtrs[row];
    }

    Logger::info("DoseInfluenceMatrix::finalize: Bucketing entries by row...");
    std::vector<size_t> rowCursor(m_rowPtrs.begin(), m_rowPtrs.begin() + m_numVoxels);
    std::vector<ColVal> rowBuckets;

    try {
        rowBuckets.resize(nnz);
    } catch (const std::bad_alloc& e) {
        Logger::error("DoseInfluenceMatrix::finalize: Cannot allocate row buckets (~" +
            std::to_string(bucketBytes / (1024ull * 1024ull * 1024ull)) +
            " GB)");
        throw std::runtime_error("Out of memory during row bucketing in sparse matrix finalization.");
    }

    for (size_t i = 0; i < nnz; ++i) {
        const size_t row = m_cooRows[i];
        const size_t pos = rowCursor[row]++;
        rowBuckets[pos] = ColVal{m_cooCols[i], m_cooVals[i]};
    }

    m_cooRows.clear();
    m_cooRows.shrink_to_fit();
    m_cooCols.clear();
    m_cooCols.shrink_to_fit();
    m_cooVals.clear();
    m_cooVals.shrink_to_fit();

    Logger::info("DoseInfluenceMatrix::finalize: Sorting columns within each row...");
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256)
#endif
    for (size_t row = 0; row < m_numVoxels; ++row) {
        const size_t begin = m_rowPtrs[row];
        const size_t end = m_rowPtrs[row + 1];
        if (end - begin <= 1) {
            continue;
        }

        std::sort(rowBuckets.begin() + static_cast<ptrdiff_t>(begin),
                  rowBuckets.begin() + static_cast<ptrdiff_t>(end),
                  [](const ColVal& lhs, const ColVal& rhs) {
                      return lhs.first < rhs.first;
                  });
    }

    Logger::info("DoseInfluenceMatrix::finalize: Counting unique entries per row...");
    std::vector<size_t> bucketRowPtrs = m_rowPtrs;
    m_rowPtrs.assign(m_numVoxels + 1, 0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256)
#endif
    for (size_t row = 0; row < m_numVoxels; ++row) {
        const size_t begin = bucketRowPtrs[row];
        const size_t end = bucketRowPtrs[row + 1];
        if (begin == end) {
            m_rowPtrs[row + 1] = 0;
            continue;
        }

        size_t uniqueCount = 1;
        for (size_t pos = begin + 1; pos < end; ++pos) {
            if (rowBuckets[pos].first != rowBuckets[pos - 1].first) {
                ++uniqueCount;
            }
        }
        m_rowPtrs[row + 1] = uniqueCount;
    }

    for (size_t row = 0; row < m_numVoxels; ++row) {
        m_rowPtrs[row + 1] += m_rowPtrs[row];
    }

    const size_t totalNnz = m_rowPtrs[m_numVoxels];
    const size_t csrBytes = totalNnz * (sizeof(size_t) + sizeof(double));
    Logger::info("DoseInfluenceMatrix::finalize: After deduplication: " +
        std::to_string(totalNnz) + " unique entries (~" +
        std::to_string(csrBytes / (1024ull * 1024ull * 1024ull)) + " GB)");

    try {
        m_values.resize(totalNnz);
        m_colIndices.resize(totalNnz);
    } catch (const std::bad_alloc& e) {
        Logger::error("DoseInfluenceMatrix::finalize: Memory allocation failed for " +
            std::to_string(totalNnz) + " entries (~" +
            std::to_string(csrBytes / (1024ull * 1024ull * 1024ull)) + " GB)");
        throw std::runtime_error("Out of memory allocating CSR storage.");
    }

    Logger::info("DoseInfluenceMatrix::finalize: Emitting CSR rows...");
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256)
#endif
    for (size_t row = 0; row < m_numVoxels; ++row) {
        const size_t begin = bucketRowPtrs[row];
        const size_t end = bucketRowPtrs[row + 1];
        size_t dst = m_rowPtrs[row];
        if (begin == end) {
            continue;
        }

        size_t currentCol = rowBuckets[begin].first;
        double currentVal = rowBuckets[begin].second;
        for (size_t pos = begin + 1; pos < end; ++pos) {
            if (rowBuckets[pos].first == currentCol) {
                currentVal += rowBuckets[pos].second;
                continue;
            }

            m_colIndices[dst] = currentCol;
            m_values[dst] = currentVal;
            ++dst;

            currentCol = rowBuckets[pos].first;
            currentVal = rowBuckets[pos].second;
        }

        m_colIndices[dst] = currentCol;
        m_values[dst] = currentVal;
    }

    Logger::info("DoseInfluenceMatrix::finalize: Freeing row bucket memory...");
    rowBuckets.clear();
    rowBuckets.shrink_to_fit();

    Logger::info("DoseInfluenceMatrix::finalize: CSR conversion complete (nnz=" +
        std::to_string(totalNnz) + ")");
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
