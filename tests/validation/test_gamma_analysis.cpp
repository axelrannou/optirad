#include <gtest/gtest.h>
#include "validation/GammaAnalysis.hpp"
#include <cmath>

namespace optirad::tests {

namespace {

// Gaussian blob centred at (cx, cy, cz) mm; peak 100 Gy.
DoseMatrix makeBlob(size_t n, double spacing, double cx, double cy, double cz,
                    double scale = 1.0, double sigma = 10.0) {
    Grid g;
    g.setDimensions(n, n, n);
    g.setSpacing(spacing, spacing, spacing);
    g.setOrigin(Vec3{0, 0, 0});
    DoseMatrix d;
    d.setGrid(g);
    d.allocate();
    for (size_t k = 0; k < n; ++k)
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0; i < n; ++i) {
                const double x = i * spacing - cx, y = j * spacing - cy, z = k * spacing - cz;
                d.at(i, j, k) = scale * 100.0 * std::exp(-(x * x + y * y + z * z) / (2 * sigma * sigma));
            }
    return d;
}

} // namespace

TEST(GammaAnalysisTest, IdenticalDosesPassCompletely) {
    auto ref = makeBlob(21, 2.0, 20, 20, 20);
    auto r = GammaAnalysis::compute(ref, ref);
    EXPECT_GT(r.numEvaluated, 0u);
    EXPECT_DOUBLE_EQ(r.passRate, 100.0);
    EXPECT_DOUBLE_EQ(r.maxGamma, 0.0);
}

TEST(GammaAnalysisTest, UniformScalingMatchesDoseCriterion) {
    auto ref = makeBlob(21, 2.0, 20, 20, 20);
    // 2% global difference at the peak is inside a 3% criterion; gamma <= ~0.67 everywhere.
    auto passing = makeBlob(21, 2.0, 20, 20, 20, 1.02);
    auto r = GammaAnalysis::compute(ref, passing);
    EXPECT_DOUBLE_EQ(r.passRate, 100.0);
    EXPECT_LE(r.maxGamma, 2.0 / 3.0 + 1e-6);

    // 6% is above the criterion at the high-dose region, with no spatial shift available to compensate.
    auto failing = makeBlob(21, 2.0, 20, 20, 20, 1.06);
    auto f = GammaAnalysis::compute(ref, failing);
    EXPECT_LT(f.passRate, 100.0);
}

TEST(GammaAnalysisTest, SmallShiftWithinDtaPasses) {
    auto ref  = makeBlob(31, 2.0, 30, 30, 30);
    auto eval = makeBlob(31, 2.0, 31, 30, 30);  // 1 mm shift, DTA = 2 mm
    auto r = GammaAnalysis::compute(ref, eval);
    EXPECT_DOUBLE_EQ(r.passRate, 100.0);
}

TEST(GammaAnalysisTest, LargeShiftFails) {
    auto ref  = makeBlob(31, 2.0, 30, 30, 30);
    auto eval = makeBlob(31, 2.0, 40, 30, 30);  // 10 mm shift
    auto r = GammaAnalysis::compute(ref, eval);
    EXPECT_LT(r.passRate, 90.0);
}

TEST(GammaAnalysisTest, DifferentGridResolutionsAreSupported) {
    auto ref  = makeBlob(31, 2.0, 30, 30, 30);
    auto eval = makeBlob(61, 1.0, 30, 30, 30);  // same physical function on a finer grid
    auto r = GammaAnalysis::compute(ref, eval);
    EXPECT_GT(r.numEvaluated, 0u);
    EXPECT_GT(r.passRate, 99.0);
}

TEST(GammaAnalysisTest, LowDoseCutoffExcludesVoxels) {
    auto ref = makeBlob(21, 2.0, 20, 20, 20);
    GammaOptions loose;  loose.lowDoseCutoffPercent = 10.0;
    GammaOptions strict; strict.lowDoseCutoffPercent = 50.0;
    EXPECT_GT(GammaAnalysis::compute(ref, ref, loose).numEvaluated,
              GammaAnalysis::compute(ref, ref, strict).numEvaluated);
}

TEST(GammaAnalysisTest, InvalidCriteriaThrow) {
    auto ref = makeBlob(5, 2.0, 4, 4, 4);
    GammaOptions bad; bad.distanceToAgreementMm = 0.0;
    EXPECT_THROW(GammaAnalysis::compute(ref, ref, bad), std::invalid_argument);
}

} // namespace optirad::tests
