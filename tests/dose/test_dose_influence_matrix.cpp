#include <gtest/gtest.h>
#include "dose/DoseInfluenceMatrix.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

namespace optirad::tests {

// ============================================================================
// Construction & basic properties
// ============================================================================

TEST(DoseInfluenceMatrixTest, DefaultConstructionEmpty) {
    DoseInfluenceMatrix dij;
    EXPECT_EQ(dij.getNumVoxels(), 0u);
    EXPECT_EQ(dij.getNumBixels(), 0u);
    EXPECT_EQ(dij.getNumNonZeros(), 0u);
    EXPECT_FALSE(dij.isFinalized());
}

TEST(DoseInfluenceMatrixTest, ConstructionSetsSize) {
    DoseInfluenceMatrix dij(100, 50);
    EXPECT_EQ(dij.getNumVoxels(), 100u);
    EXPECT_EQ(dij.getNumBixels(), 50u);
}

TEST(DoseInfluenceMatrixTest, SetDimensionsClearsState) {
    DoseInfluenceMatrix dij(10, 5);
    dij.setValue(0, 0, 1.0);
    EXPECT_EQ(dij.getNumNonZeros(), 1u);

    dij.setDimensions(20, 10);
    EXPECT_EQ(dij.getNumVoxels(), 20u);
    EXPECT_EQ(dij.getNumBixels(), 10u);
    EXPECT_EQ(dij.getNumNonZeros(), 0u);
    EXPECT_FALSE(dij.isFinalized());
}

// ============================================================================
// setValue / getValue (COO phase)
// ============================================================================

TEST(DoseInfluenceMatrixTest, SetValueBeforeFinalize) {
    DoseInfluenceMatrix dij(10, 5);
    dij.setValue(3, 2, 4.5);
    EXPECT_EQ(dij.getNumNonZeros(), 1u);

    // getValue before finalize (linear scan)
    EXPECT_DOUBLE_EQ(dij.getValue(3, 2), 4.5);
    EXPECT_DOUBLE_EQ(dij.getValue(0, 0), 0.0);
}

TEST(DoseInfluenceMatrixTest, SetValueOutOfBoundsThrows) {
    DoseInfluenceMatrix dij(10, 5);
    EXPECT_THROW(dij.setValue(10, 0, 1.0), std::out_of_range);
    EXPECT_THROW(dij.setValue(0, 5, 1.0), std::out_of_range);
}

TEST(DoseInfluenceMatrixTest, SetValueAfterFinalizeThrows) {
    DoseInfluenceMatrix dij(10, 5);
    dij.setValue(0, 0, 1.0);
    dij.finalize();
    EXPECT_THROW(dij.setValue(1, 1, 2.0), std::runtime_error);
}

TEST(DoseInfluenceMatrixTest, DuplicateEntriesSummed) {
    DoseInfluenceMatrix dij(5, 3);
    dij.setValue(1, 2, 3.0);
    dij.setValue(1, 2, 7.0); // duplicate → should sum to 10.0
    dij.finalize();
    EXPECT_DOUBLE_EQ(dij.getValue(1, 2), 10.0);
}

// ============================================================================
// Finalize (COO → CSR)
// ============================================================================

TEST(DoseInfluenceMatrixTest, FinalizeEmpty) {
    DoseInfluenceMatrix dij(10, 5);
    dij.finalize(); // no entries
    EXPECT_TRUE(dij.isFinalized());
    EXPECT_EQ(dij.getNumNonZeros(), 0u);
    EXPECT_EQ(dij.getRowPtrs().size(), 11u);
}

TEST(DoseInfluenceMatrixTest, FinalizeProducesValidCSR) {
    DoseInfluenceMatrix dij(4, 3);
    //  Row 0: (0, 1)=2.0
    //  Row 1: (1, 0)=3.0, (1, 2)=4.0
    //  Row 2: empty
    //  Row 3: (3, 1)=5.0
    dij.setValue(0, 1, 2.0);
    dij.setValue(1, 0, 3.0);
    dij.setValue(1, 2, 4.0);
    dij.setValue(3, 1, 5.0);
    dij.finalize();

    EXPECT_TRUE(dij.isFinalized());
    EXPECT_EQ(dij.getNumNonZeros(), 4u);

    const auto& rowPtrs = dij.getRowPtrs();
    EXPECT_EQ(rowPtrs.size(), 5u);
    EXPECT_EQ(rowPtrs[0], 0u); // row 0 starts at 0
    EXPECT_EQ(rowPtrs[1], 1u); // row 1 starts at 1
    EXPECT_EQ(rowPtrs[2], 3u); // row 2 starts at 3
    EXPECT_EQ(rowPtrs[3], 3u); // row 3 starts at 3 (row 2 empty)
    EXPECT_EQ(rowPtrs[4], 4u); // end

    // Verify values
    EXPECT_DOUBLE_EQ(dij.getValue(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(dij.getValue(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(dij.getValue(1, 2), 4.0);
    EXPECT_DOUBLE_EQ(dij.getValue(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(dij.getValue(3, 1), 5.0);
}

TEST(DoseInfluenceMatrixTest, DoubleFinalizeSafe) {
    DoseInfluenceMatrix dij(5, 3);
    dij.setValue(0, 0, 1.0);
    dij.finalize();
    dij.finalize(); // should be a no-op
    EXPECT_TRUE(dij.isFinalized());
    EXPECT_DOUBLE_EQ(dij.getValue(0, 0), 1.0);
}

// ============================================================================
// computeDose
// ============================================================================

TEST(DoseInfluenceMatrixTest, ComputeDoseIdentity) {
    // 3 voxels, 3 bixels, diagonal Dij → dose = weights
    DoseInfluenceMatrix dij(3, 3);
    dij.setValue(0, 0, 1.0);
    dij.setValue(1, 1, 1.0);
    dij.setValue(2, 2, 1.0);
    dij.finalize();

    std::vector<double> weights = {10.0, 20.0, 30.0};
    auto dose = dij.computeDose(weights);
    ASSERT_EQ(dose.size(), 3u);
    EXPECT_DOUBLE_EQ(dose[0], 10.0);
    EXPECT_DOUBLE_EQ(dose[1], 20.0);
    EXPECT_DOUBLE_EQ(dose[2], 30.0);
}

TEST(DoseInfluenceMatrixTest, ComputeDoseSumMultipleBixels) {
    // 2 voxels, 3 bixels
    DoseInfluenceMatrix dij(2, 3);
    dij.setValue(0, 0, 1.0);
    dij.setValue(0, 1, 2.0);
    dij.setValue(0, 2, 3.0);
    dij.setValue(1, 1, 4.0);
    dij.finalize();

    std::vector<double> weights = {1.0, 1.0, 1.0};
    auto dose = dij.computeDose(weights);
    ASSERT_EQ(dose.size(), 2u);
    EXPECT_DOUBLE_EQ(dose[0], 6.0); // 1+2+3
    EXPECT_DOUBLE_EQ(dose[1], 4.0); // 0+4+0
}

TEST(DoseInfluenceMatrixTest, ComputeDoseWithZeroWeights) {
    DoseInfluenceMatrix dij(2, 2);
    dij.setValue(0, 0, 5.0);
    dij.setValue(1, 1, 3.0);
    dij.finalize();

    std::vector<double> weights = {0.0, 0.0};
    auto dose = dij.computeDose(weights);
    EXPECT_DOUBLE_EQ(dose[0], 0.0);
    EXPECT_DOUBLE_EQ(dose[1], 0.0);
}

TEST(DoseInfluenceMatrixTest, ComputeDoseBeforeFinalizeThrows) {
    DoseInfluenceMatrix dij(3, 3);
    dij.setValue(0, 0, 1.0);
    std::vector<double> weights = {1.0, 1.0, 1.0};
    EXPECT_THROW(dij.computeDose(weights), std::runtime_error);
}

// ============================================================================
// accumulateTransposeProduct
// ============================================================================

TEST(DoseInfluenceMatrixTest, TransposeProductDiagonal) {
    DoseInfluenceMatrix dij(3, 3);
    dij.setValue(0, 0, 2.0);
    dij.setValue(1, 1, 3.0);
    dij.setValue(2, 2, 4.0);
    dij.finalize();

    std::vector<double> voxelGrad = {1.0, 1.0, 1.0};
    std::vector<double> grad(3, 0.0);
    dij.accumulateTransposeProduct(voxelGrad, grad);

    EXPECT_DOUBLE_EQ(grad[0], 2.0);
    EXPECT_DOUBLE_EQ(grad[1], 3.0);
    EXPECT_DOUBLE_EQ(grad[2], 4.0);
}

TEST(DoseInfluenceMatrixTest, TransposeProductMultipleVoxels) {
    // 2 voxels, 2 bixels
    // dij = [[1, 2], [3, 4]]
    // voxelGrad = [1, 1]
    // grad += dij^T * voxelGrad = [1+3, 2+4] = [4, 6]
    DoseInfluenceMatrix dij(2, 2);
    dij.setValue(0, 0, 1.0);
    dij.setValue(0, 1, 2.0);
    dij.setValue(1, 0, 3.0);
    dij.setValue(1, 1, 4.0);
    dij.finalize();

    std::vector<double> voxelGrad = {1.0, 1.0};
    std::vector<double> grad(2, 0.0);
    dij.accumulateTransposeProduct(voxelGrad, grad);

    EXPECT_DOUBLE_EQ(grad[0], 4.0);
    EXPECT_DOUBLE_EQ(grad[1], 6.0);
}

// ============================================================================
// loadCSR
// ============================================================================

TEST(DoseInfluenceMatrixTest, LoadCSRDirectly) {
    DoseInfluenceMatrix dij;
    dij.setDimensions(3, 4);

    // Row 0: (0,1)=2.5
    // Row 1: (1,0)=1.0 (1,3)=3.0
    // Row 2: empty
    std::vector<size_t> rowPtrs = {0, 1, 3, 3};
    std::vector<size_t> colIndices = {1, 0, 3};
    std::vector<double> values = {2.5, 1.0, 3.0};

    dij.loadCSR(std::move(rowPtrs), std::move(colIndices), std::move(values));

    EXPECT_TRUE(dij.isFinalized());
    EXPECT_EQ(dij.getNumNonZeros(), 3u);
    EXPECT_DOUBLE_EQ(dij.getValue(0, 1), 2.5);
    EXPECT_DOUBLE_EQ(dij.getValue(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(dij.getValue(1, 3), 3.0);
    EXPECT_DOUBLE_EQ(dij.getValue(2, 0), 0.0);
}

// ============================================================================
// Large sparse matrix (memory efficiency test)
// ============================================================================

TEST(DoseInfluenceMatrixTest, LargeSparseDoesNotOOM) {
    // 1 million voxels, 1000 bixels — only 5000 non-zeros
    // With the old dense approach this would need ~8 GB. Now it uses ~kilobytes.
    size_t numVoxels = 1000000;
    size_t numBixels = 1000;
    DoseInfluenceMatrix dij(numVoxels, numBixels);

    // Insert 5000 entries
    for (size_t i = 0; i < 5000; ++i) {
        dij.setValue(i * 200, i % numBixels, 0.01 * (i + 1));
    }
    dij.finalize();

    EXPECT_EQ(dij.getNumNonZeros(), 5000u);
    EXPECT_DOUBLE_EQ(dij.getValue(0, 0), 0.01);
    EXPECT_DOUBLE_EQ(dij.getValue(200, 1), 0.02);

    // Compute dose with unit weights
    std::vector<double> weights(numBixels, 1.0);
    auto dose = dij.computeDose(weights);
    EXPECT_EQ(dose.size(), numVoxels);

    // Verify one known entry
    EXPECT_DOUBLE_EQ(dose[0], 0.01);
}

// ============================================================================
// reserveNonZeros
// ============================================================================

TEST(DoseInfluenceMatrixTest, ReserveDoesNotAffectResult) {
    DoseInfluenceMatrix dij(5, 5);
    dij.reserveNonZeros(100); // hint, should not change behavior
    dij.setValue(0, 0, 1.0);
    dij.setValue(4, 4, 2.0);
    dij.finalize();
    EXPECT_DOUBLE_EQ(dij.getValue(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(dij.getValue(4, 4), 2.0);
}

// ============================================================================
// operator() const (alias for getValue)
// ============================================================================

TEST(DoseInfluenceMatrixTest, ConstOperatorWorks) {
    DoseInfluenceMatrix dij(5, 5);
    dij.setValue(2, 3, 7.7);
    dij.finalize();

    const DoseInfluenceMatrix& cdij = dij;
    EXPECT_DOUBLE_EQ(cdij(2, 3), 7.7);
    EXPECT_DOUBLE_EQ(cdij(0, 0), 0.0);
}

} // namespace optirad::tests
