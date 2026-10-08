#include "DoseInfluenceMatrix.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace optirad {

namespace {

constexpr size_t kMaxIndex = std::numeric_limits<uint32_t>::max();

void checkDimensions(size_t numVoxels, size_t numBixels) {
    if (numVoxels >= kMaxIndex || numBixels >= kMaxIndex)
        throw std::length_error("DoseInfluenceMatrix: dimensions must fit in uint32");
}

} // namespace

// ────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────

DoseInfluenceMatrix::DoseInfluenceMatrix(size_t numVoxels, size_t numBixels)
    : m_numVoxels(numVoxels), m_numBixels(numBixels) {
    checkDimensions(numVoxels, numBixels);
}

void DoseInfluenceMatrix::setDimensions(size_t numVoxels, size_t numBixels) {
    checkDimensions(numVoxels, numBixels);
    m_numVoxels = numVoxels;
    m_numBixels = numBixels;
    m_finalized = false;
    m_cooRows.clear();
    m_cooCols.clear();
    m_cooVals.clear();
    m_blocks.clear();
    m_rowMask.clear();
}

void DoseInfluenceMatrix::setComputedRows(std::vector<uint8_t> mask) {
    if (!mask.empty() && mask.size() != m_numVoxels)
        throw std::invalid_argument("DoseInfluenceMatrix::setComputedRows: mask size mismatch");
    m_rowMask = std::move(mask);
}

void DoseInfluenceMatrix::reserveNonZeros(size_t nnz) {
    m_cooRows.reserve(nnz);
    m_cooCols.reserve(nnz);
    m_cooVals.reserve(nnz);
}

void DoseInfluenceMatrix::appendBatch(const std::vector<uint32_t>& rows,
                                      const std::vector<uint32_t>& cols,
                                      const std::vector<float>& vals) {
    if (m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix: cannot appendBatch after finalize()");
    m_cooRows.insert(m_cooRows.end(), rows.begin(), rows.end());
    m_cooCols.insert(m_cooCols.end(), cols.begin(), cols.end());
    m_cooVals.insert(m_cooVals.end(), vals.begin(), vals.end());
}

void DoseInfluenceMatrix::setValue(size_t voxel, size_t bixel, double value) {
    if (m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix: cannot setValue after finalize()");
    if (voxel >= m_numVoxels || bixel >= m_numBixels)
        throw std::out_of_range("DoseInfluenceMatrix::setValue: index out of bounds");
    m_cooRows.push_back(static_cast<uint32_t>(voxel));
    m_cooCols.push_back(static_cast<uint32_t>(bixel));
    m_cooVals.push_back(static_cast<float>(value));
}

// ────────────────────────────────────────────────────────────────
// COO → CSR block
// ────────────────────────────────────────────────────────────────

void DoseInfluenceMatrix::endBlock() {
    if (m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix: cannot endBlock after finalize()");

    const size_t nnz = m_cooRows.size();
    if (nnz == 0) return;
    if (nnz >= kMaxIndex)
        throw std::length_error("DoseInfluenceMatrix::endBlock: block has too many entries (>= 2^32)");

    using ColVal = std::pair<uint32_t, float>;

    std::vector<uint32_t> bucketPtrs(m_numVoxels + 1, 0);
    for (size_t i = 0; i < nnz; ++i) {
        if (m_cooRows[i] >= m_numVoxels || m_cooCols[i] >= m_numBixels)
            throw std::runtime_error("DoseInfluenceMatrix::endBlock: index out of bounds at entry " +
                                     std::to_string(i));
        ++bucketPtrs[m_cooRows[i] + 1];
    }
    for (size_t row = 0; row < m_numVoxels; ++row) bucketPtrs[row + 1] += bucketPtrs[row];

    std::vector<ColVal> buckets(nnz);
    {
        std::vector<uint32_t> cursor(bucketPtrs.begin(), bucketPtrs.end() - 1);
        for (size_t i = 0; i < nnz; ++i)
            buckets[cursor[m_cooRows[i]]++] = ColVal{m_cooCols[i], m_cooVals[i]};
    }
    std::vector<uint32_t>().swap(m_cooRows);
    std::vector<uint32_t>().swap(m_cooCols);
    std::vector<float>().swap(m_cooVals);

    Block block;
    block.rowPtrs.assign(m_numVoxels + 1, 0);

    // Sort each row by column and count unique columns.
    #pragma omp parallel for schedule(dynamic, 256)
    for (size_t row = 0; row < m_numVoxels; ++row) {
        const size_t begin = bucketPtrs[row], end = bucketPtrs[row + 1];
        if (begin == end) continue;
        std::sort(buckets.begin() + static_cast<ptrdiff_t>(begin),
                  buckets.begin() + static_cast<ptrdiff_t>(end),
                  [](const ColVal& a, const ColVal& b) { return a.first < b.first; });
        uint32_t unique = 1;
        for (size_t p = begin + 1; p < end; ++p)
            if (buckets[p].first != buckets[p - 1].first) ++unique;
        block.rowPtrs[row + 1] = unique;
    }
    for (size_t row = 0; row < m_numVoxels; ++row) block.rowPtrs[row + 1] += block.rowPtrs[row];

    const size_t total = block.rowPtrs[m_numVoxels];
    block.colIndices.resize(total);
    block.values.resize(total);

    // Emit rows, summing duplicate columns.
    #pragma omp parallel for schedule(dynamic, 256)
    for (size_t row = 0; row < m_numVoxels; ++row) {
        const size_t begin = bucketPtrs[row], end = bucketPtrs[row + 1];
        if (begin == end) continue;
        size_t dst = block.rowPtrs[row];
        uint32_t col = buckets[begin].first;
        double sum = buckets[begin].second;
        for (size_t p = begin + 1; p < end; ++p) {
            if (buckets[p].first == col) { sum += buckets[p].second; continue; }
            block.colIndices[dst] = col;
            block.values[dst++] = static_cast<float>(sum);
            col = buckets[p].first;
            sum = buckets[p].second;
        }
        block.colIndices[dst] = col;
        block.values[dst] = static_cast<float>(sum);
    }

    m_blocks.push_back(std::move(block));
}

void DoseInfluenceMatrix::finalize() {
    if (m_finalized) return;
    endBlock();
    if (m_blocks.empty()) {
        Logger::warn("DoseInfluenceMatrix::finalize: No entries to finalize");
        Block empty;
        empty.rowPtrs.assign(m_numVoxels + 1, 0);
        m_blocks.push_back(std::move(empty));
    }
    m_finalized = true;
}

// ────────────────────────────────────────────────────────────────
// Direct loading (deserialization)
// ────────────────────────────────────────────────────────────────

void DoseInfluenceMatrix::addBlock(Block block) {
    if (block.rowPtrs.size() != m_numVoxels + 1)
        throw std::invalid_argument("DoseInfluenceMatrix::addBlock: rowPtrs size mismatch");
    m_blocks.push_back(std::move(block));
    m_finalized = true;
}

void DoseInfluenceMatrix::loadCSR(std::vector<size_t> rowPtrs,
                                  std::vector<size_t> colIndices,
                                  std::vector<double>  values) {
    if (colIndices.size() != values.size() || colIndices.size() >= kMaxIndex)
        throw std::invalid_argument("DoseInfluenceMatrix::loadCSR: inconsistent or oversized arrays");
    Block block;
    block.rowPtrs.assign(rowPtrs.begin(), rowPtrs.end());
    block.colIndices.assign(colIndices.begin(), colIndices.end());
    block.values.assign(values.begin(), values.end());
    m_blocks.clear();
    m_cooRows.clear(); m_cooCols.clear(); m_cooVals.clear();
    addBlock(std::move(block));
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

    double sum = 0.0;
    // Pending COO entries (slow linear scan; only for debugging/tests).
    for (size_t k = 0; k < m_cooRows.size(); ++k)
        if (m_cooRows[k] == voxel && m_cooCols[k] == bixel) sum += m_cooVals[k];

    for (const auto& b : m_blocks) {
        const auto begin = b.colIndices.begin() + b.rowPtrs[voxel];
        const auto end   = b.colIndices.begin() + b.rowPtrs[voxel + 1];
        const auto it = std::lower_bound(begin, end, static_cast<uint32_t>(bixel));
        if (it != end && *it == bixel) sum += b.values[static_cast<size_t>(it - b.colIndices.begin())];
    }
    return sum;
}

size_t DoseInfluenceMatrix::getNumVoxels() const { return m_numVoxels; }
size_t DoseInfluenceMatrix::getNumBixels() const { return m_numBixels; }

size_t DoseInfluenceMatrix::getNumNonZeros() const {
    size_t n = m_cooRows.size();
    for (const auto& b : m_blocks) n += b.values.size();
    return n;
}

double DoseInfluenceMatrix::getMaxValue() const {
    if (!m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix::getMaxValue requires finalize()");
    float best = 0.0f;
    for (const auto& b : m_blocks)
        if (!b.values.empty()) best = std::max(best, *std::max_element(b.values.begin(), b.values.end()));
    return best;
}

// ────────────────────────────────────────────────────────────────
// Linear algebra (require finalized)
// ────────────────────────────────────────────────────────────────

std::vector<double> DoseInfluenceMatrix::computeDose(const std::vector<double>& weights) const {
    if (!m_finalized)
        throw std::runtime_error("DoseInfluenceMatrix::computeDose requires finalize()");

    std::vector<double> dose(m_numVoxels, 0.0);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t v = 0; v < m_numVoxels; ++v) {
        double sum = 0.0;
        for (const auto& b : m_blocks)
            for (size_t k = b.rowPtrs[v]; k < b.rowPtrs[v + 1]; ++k)
                sum += static_cast<double>(b.values[k]) * weights[b.colIndices[k]];
        dose[v] = sum;
    }
    return dose;
}

void DoseInfluenceMatrix::accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                                     std::vector<double>& grad) const {
    if (!m_finalized)
        throw std::runtime_error("accumulateTransposeProduct requires finalize()");

    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t v = 0; v < m_numVoxels; ++v) {
        const double gv = voxelGrad[v];
        if (gv == 0.0) continue;
        for (const auto& b : m_blocks) {
            for (size_t k = b.rowPtrs[v]; k < b.rowPtrs[v + 1]; ++k) {
                #pragma omp atomic
                grad[b.colIndices[k]] += gv * static_cast<double>(b.values[k]);
            }
        }
    }
}

} // namespace optirad
