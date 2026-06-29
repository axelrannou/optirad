#include <gtest/gtest.h>
#include "sequencing/MlcLeafBounds.hpp"
#include "core/Machine.hpp"

using namespace optirad;

// ─────────────────────── Helpers ───────────────────────

static MachineGeometry makeMlcGeometry(
    int numLeaves,
    const std::vector<double>& leafWidths,
    int numInnerPairs)
{
    MachineGeometry g;
    g.numLeaves     = numLeaves;
    g.leafWidths    = leafWidths;
    g.numInnerPairs = numInnerPairs;
    return g;
}

static MachineGeometry millennium120()
{
    // Varian Millennium 120 MLC
    // 120 leaves / 60 pairs: 10 outer (10mm) + 40 inner (5mm) + 10 outer (10mm)
    // Total span: 10*10 + 40*5 + 10*10 = 400mm  (±200mm)
    return makeMlcGeometry(120, {5.0, 10.0}, 40);
}

static MachineGeometry hd120()
{
    // Varian HD120 MLC
    // 120 leaves / 60 pairs: 14 outer (5mm) + 32 inner (2.5mm) + 14 outer (5mm)
    // Total span: 14*5 + 32*2.5 + 14*5 = 220mm  (±110mm)
    return makeMlcGeometry(120, {2.5, 5.0}, 32);
}

// ═══════════════════════════════════════════════════════════════
//  Millennium 120  tests
// ═══════════════════════════════════════════════════════════════

TEST(Millennium120, BoundaryCount)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    // 60 pairs → 61 boundaries
    ASSERT_EQ(static_cast<int>(bounds.size()), 61);
}

TEST(Millennium120, TotalSpan)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    EXPECT_DOUBLE_EQ(bounds.front(), -200.0);
    EXPECT_DOUBLE_EQ(bounds.back(),  +200.0);
}

TEST(Millennium120, Symmetry)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    // Centered at 0: bounds[0] + bounds[60] should be 0
    EXPECT_DOUBLE_EQ(bounds[0] + bounds[60], 0.0);
}

TEST(Millennium120, OuterPairWidth)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    // Pairs 0–9 are outer (10mm each)
    for (int lp = 0; lp < 10; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 10.0)
            << "outer pair " << lp;
    }
}

TEST(Millennium120, InnerPairWidth)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    // Pairs 10–49 are inner (5mm each)
    for (int lp = 10; lp < 50; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 5.0)
            << "inner pair " << lp;
    }
}

TEST(Millennium120, TopOuterPairWidth)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    // Pairs 50–59 are outer (10mm each)
    for (int lp = 50; lp < 60; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 10.0)
            << "top-outer pair " << lp;
    }
}

TEST(Millennium120, SpecificBoundaries)
{
    auto bounds = buildMlcLeafBounds(millennium120());
    //  bounds[0]  = -200.0  (bottom edge of pair 0)
    //  bounds[10] = -100.0  (transition outer→inner)
    //  bounds[30] =    0.0  (midpoint = center of pair 29/30)
    //  bounds[50] = +100.0  (transition inner→outer)
    //  bounds[60] = +200.0  (top edge of pair 59)
    EXPECT_DOUBLE_EQ(bounds[0],   -200.0);
    EXPECT_DOUBLE_EQ(bounds[10],  -100.0);
    EXPECT_DOUBLE_EQ(bounds[30],     0.0);
    EXPECT_DOUBLE_EQ(bounds[50],  +100.0);
    EXPECT_DOUBLE_EQ(bounds[60],  +200.0);
}

// ═══════════════════════════════════════════════════════════════
//  HD120  tests
// ═══════════════════════════════════════════════════════════════

TEST(HD120, BoundaryCount)
{
    auto bounds = buildMlcLeafBounds(hd120());
    // 60 pairs → 61 boundaries
    ASSERT_EQ(static_cast<int>(bounds.size()), 61);
}

TEST(HD120, TotalSpan)
{
    auto bounds = buildMlcLeafBounds(hd120());
    // 14*5 + 32*2.5 + 14*5 = 70 + 80 + 70 = 220mm → ±110mm
    EXPECT_DOUBLE_EQ(bounds.front(), -110.0);
    EXPECT_DOUBLE_EQ(bounds.back(),  +110.0);
}

TEST(HD120, Symmetry)
{
    auto bounds = buildMlcLeafBounds(hd120());
    EXPECT_DOUBLE_EQ(bounds[0] + bounds[60], 0.0);
}

TEST(HD120, OuterPairWidth)
{
    auto bounds = buildMlcLeafBounds(hd120());
    // Pairs 0–13 are outer (5mm each)
    for (int lp = 0; lp < 14; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 5.0)
            << "outer pair " << lp;
    }
}

TEST(HD120, InnerPairWidth)
{
    auto bounds = buildMlcLeafBounds(hd120());
    // Pairs 14–45 are inner (2.5mm each)
    for (int lp = 14; lp < 46; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 2.5)
            << "inner pair " << lp;
    }
}

TEST(HD120, TopOuterPairWidth)
{
    auto bounds = buildMlcLeafBounds(hd120());
    // Pairs 46–59 are outer (5mm each)
    for (int lp = 46; lp < 60; ++lp) {
        EXPECT_DOUBLE_EQ(bounds[lp + 1] - bounds[lp], 5.0)
            << "top-outer pair " << lp;
    }
}

TEST(HD120, SpecificBoundaries)
{
    auto bounds = buildMlcLeafBounds(hd120());
    //  bounds[0]  = -110.0  (bottom edge of pair 0)
    //  bounds[14] =  -30.0  (14*5 - 110 = 70-110 = -40? no: -110+14*5 = -110+70 = -40)
    //  bounds[30] =    0.0  midpoint check (14*5 + 16*2.5 - 110 = 70+40-110 = 0)
    //  bounds[46] =  +40.0  (transition inner→outer: -110 + 14*5 + 32*2.5 = -110+70+80 = +40)
    //  bounds[60] = +110.0
    EXPECT_DOUBLE_EQ(bounds[0],   -110.0);
    EXPECT_DOUBLE_EQ(bounds[14],   -40.0);
    EXPECT_DOUBLE_EQ(bounds[30],     0.0);
    EXPECT_DOUBLE_EQ(bounds[46],   +40.0);
    EXPECT_DOUBLE_EQ(bounds[60],  +110.0);
}

// ═══════════════════════════════════════════════════════════════
//  Machine switching tests
// ═══════════════════════════════════════════════════════════════

TEST(MachineSwitching, DifferentSpans)
{
    auto m120 = buildMlcLeafBounds(millennium120());
    auto hd   = buildMlcLeafBounds(hd120());

    EXPECT_NE(m120.front(), hd.front());
    EXPECT_NE(m120.back(),  hd.back());
    EXPECT_DOUBLE_EQ(m120.front(), -200.0);
    EXPECT_DOUBLE_EQ(hd.front(),   -110.0);
}

TEST(MachineSwitching, SamePairCount)
{
    auto m120 = buildMlcLeafBounds(millennium120());
    auto hd   = buildMlcLeafBounds(hd120());
    // Both have 120 leaves → 60 pairs → 61 boundaries
    EXPECT_EQ(m120.size(), hd.size());
}

TEST(MachineSwitching, SwitchChangesInnerWidth)
{
    auto m120 = buildMlcLeafBounds(millennium120());
    auto hd   = buildMlcLeafBounds(hd120());

    double m120InnerWidth = m120[11] - m120[10];  // inner pair
    double hdInnerWidth   = hd[15]   - hd[14];    // inner pair

    EXPECT_DOUBLE_EQ(m120InnerWidth, 5.0);
    EXPECT_DOUBLE_EQ(hdInnerWidth,   2.5);
}

// ═══════════════════════════════════════════════════════════════
//  Edge cases
// ═══════════════════════════════════════════════════════════════

TEST(MlcLeafBounds, EmptyGeometry)
{
    MachineGeometry g;  // default: numLeaves=0
    auto bounds = buildMlcLeafBounds(g);
    EXPECT_TRUE(bounds.empty());
}

TEST(MlcLeafBounds, UniformWidth)
{
    // Uniform-width MLC: 10 leaves, 5mm each → 5 pairs
    MachineGeometry g = makeMlcGeometry(10, {5.0}, 0);
    auto bounds = buildMlcLeafBounds(g);
    ASSERT_EQ(static_cast<int>(bounds.size()), 6);
    EXPECT_DOUBLE_EQ(bounds.front(), -12.5);
    EXPECT_DOUBLE_EQ(bounds.back(),  +12.5);
    for (int i = 0; i < 5; ++i)
        EXPECT_DOUBLE_EQ(bounds[i + 1] - bounds[i], 5.0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
