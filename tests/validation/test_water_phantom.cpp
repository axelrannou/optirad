#include <gtest/gtest.h>
#include "validation/WaterPhantom.hpp"
#include <algorithm>

namespace optirad::tests {

namespace {

// Computed once for all tests (about 15 s). No Dij thresholding, as in matRad.
class WaterPhantomBenchmark : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        WaterPhantomConfig cfg;
        cfg.relativeThreshold = 0.0;
        result = new WaterPhantomResult(WaterPhantom::compute(cfg));
    }
    static void TearDownTestSuite() { delete result; result = nullptr; }

    static WaterPhantomResult* result;
};
WaterPhantomResult* WaterPhantomBenchmark::result = nullptr;

// Positive x position where the profile at height y falls through 50 % of its central value.
double rightEdge(const DoseMatrix& dose, double y) {
    const auto p = sampleProfileX(dose, y, 0, 0, 120, 0.25);
    const double half = 0.5 * p.dose.front();
    for (size_t i = 0; i < p.dose.size(); ++i)
        if (p.dose[i] < half) return p.position[i];
    return p.position.back();
}

} // namespace

// Published 6 MV, 10x10 cm2 values (BJR Suppl. 25; TRS-398 tables): PDD(10 cm) ~ 66-67 %, PDD(20 cm) ~ 38-39 %
// at SSD 100 cm. The phantom uses SSD 90 cm, which lowers the PDDs by about 1 point; tolerances are generous.
TEST_F(WaterPhantomBenchmark, DepthDoseMatchesPublished6MVRange) {
    const auto c = sampleDepthDose(*result->dose, 0, 0, result->surfaceY, result->surfaceY + 300.0, 1.0);
    const auto it = std::max_element(c.dose.begin(), c.dose.end());
    const double dmaxDepth = c.position[it - c.dose.begin()] - result->surfaceY;

    EXPECT_GE(dmaxDepth, 10.0);
    EXPECT_LE(dmaxDepth, 20.0);
    EXPECT_NEAR(100.0 * c.dose[100] / *it, 65.0, 3.0);
    EXPECT_NEAR(100.0 * c.dose[200] / *it, 37.0, 3.0);
}

// The bixel grid extends the 10 cm target by the STF margin, so the 50 % edge sits slightly outside +-50 mm.
TEST_F(WaterPhantomBenchmark, ProfileEdgesAreSymmetricAndNearFieldBorder) {
    const auto p = sampleProfileX(*result->dose, 0.0, 0, -100, 100, 0.5);
    const double central = p.dose[p.dose.size() / 2];

    double left = 0.0, right = 0.0;
    for (size_t i = 0; i < p.dose.size(); ++i)
        if (p.dose[i] >= 0.5 * central) { left = p.position[i]; break; }
    for (size_t i = p.dose.size(); i-- > 0;)
        if (p.dose[i] >= 0.5 * central) { right = p.position[i]; break; }

    EXPECT_NEAR(left + right, 0.0, 3.0);
    EXPECT_GT(right, 48.0);
    EXPECT_LT(right, 65.0);
}

// Regression: off-axis bixels once stayed parallel to the beam axis instead of diverging from the source.
TEST_F(WaterPhantomBenchmark, FieldDivergesWithDepth) {
    const double shallow = rightEdge(*result->dose, -50.0);  // 50 mm deep, beam magnification (1000-50)/1000
    const double deep    = rightEdge(*result->dose, 100.0);  // 200 mm deep, magnification 1.10
    EXPECT_NEAR(deep / shallow, 1.10 / 0.95, 0.05);
}

} // namespace optirad::tests

namespace optirad::tests {

namespace {

WaterPhantomConfig smallConfig() {
    WaterPhantomConfig c;
    c.resolutionMm = 5.0;
    c.bixelWidthMm = 5.0;
    c.fieldSizeMm = 40.0;
    c.lateralHalfSizeMm = 60.0;
    c.depthMm = 150.0;
    c.relativeThreshold = 0.0;
    c.computeDirect = true;
    return c;
}

} // namespace

// The direct dose must equal Dij x weights (up to the float storage of the Dij).
TEST(DirectDoseTest, MatchesDijDose) {
    const auto r = WaterPhantom::compute(smallConfig());
    ASSERT_TRUE(r.dijHasAllRows);
    ASSERT_EQ(r.directDoses.size(), 1u);
    const auto& a = *r.doses[0];
    const auto& b = *r.directDoses[0];
    ASSERT_EQ(a.size(), b.size());
    const double max = a.getMax();
    ASSERT_GT(max, 0.0);
    for (size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a.data()[i], b.data()[i], 1e-5 * max);
}

// Skipping the body drops Dij rows, but the direct dose still covers the whole phantom.
TEST(DirectDoseTest, ExternalSkippedInDijButNotInDirectDose) {
    auto cfg = smallConfig();
    cfg.excludeExternal = true;
    const auto restricted = WaterPhantom::compute(cfg);
    cfg.excludeExternal = false;
    const auto full = WaterPhantom::compute(cfg);

    EXPECT_FALSE(restricted.dijHasAllRows);
    EXPECT_LT(restricted.dijNonZeros, full.dijNonZeros);

    // Direct dose does not depend on the Dij restriction.
    const auto& d = *restricted.directDoses[0];
    const auto& f = *full.doses[0];
    const double max = f.getMax();
    for (size_t i = 0; i < f.size(); ++i) EXPECT_NEAR(d.data()[i], f.data()[i], 1e-5 * max);

    // Where the restricted Dij has rows, it agrees with the full one.
    const auto& r = *restricted.doses[0];
    size_t covered = 0;
    for (size_t i = 0; i < r.size(); ++i)
        if (r.data()[i] > 0.0) { ++covered; EXPECT_NEAR(r.data()[i], f.data()[i], 1e-5 * max); }
    EXPECT_GT(covered, 0u);
    EXPECT_LT(covered, f.size());
}

} // namespace optirad::tests
