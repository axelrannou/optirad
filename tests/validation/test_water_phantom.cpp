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
