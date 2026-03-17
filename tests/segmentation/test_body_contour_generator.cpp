#include <gtest/gtest.h>
#include "segmentation/BodyContourGenerator.hpp"
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"
#include "geometry/Structure.hpp"
#include <cmath>

using namespace optirad;

namespace {

/**
 * Helper: create a small CT volume with given dimensions and fill with a default HU value.
 */
std::unique_ptr<Volume<int16_t>> createTestVolume(size_t nx, size_t ny, size_t nz,
                                                   int16_t fillHU = -1000,
                                                   double spacingX = 1.0,
                                                   double spacingY = 1.0,
                                                   double spacingZ = 1.0) {
    auto vol = std::make_unique<Volume<int16_t>>();
    Grid grid;
    grid.setDimensions(nx, ny, nz);
    grid.setSpacing(spacingX, spacingY, spacingZ);
    grid.setOrigin({0.0, 0.0, 0.0});
    vol->setGrid(grid);
    vol->allocate();

    // Fill with background HU
    int16_t* data = vol->data();
    for (size_t i = 0; i < vol->size(); ++i) {
        data[i] = fillHU;
    }
    return vol;
}

/**
 * Helper: place a solid sphere of given HU into a volume.
 */
void placeSphere(Volume<int16_t>& vol, double cx, double cy, double cz,
                 double radius, int16_t huValue) {
    auto dims = vol.getGrid().getDimensions();
    size_t nx = dims[0], ny = dims[1], nz = dims[2];
    int16_t* data = vol.data();

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                double dx = static_cast<double>(i) - cx;
                double dy = static_cast<double>(j) - cy;
                double dz = static_cast<double>(k) - cz;
                if (dx*dx + dy*dy + dz*dz <= radius*radius) {
                    data[k * ny * nx + j * nx + i] = huValue;
                }
            }
        }
    }
}

/**
 * Helper: place a hollow cylinder (outer shell only) into a volume.
 * Oriented along Z axis.
 */
void placeHollowCylinder(Volume<int16_t>& vol, double cx, double cy,
                          double outerRadius, double innerRadius,
                          int16_t huValue) {
    auto dims = vol.getGrid().getDimensions();
    size_t nx = dims[0], ny = dims[1], nz = dims[2];
    int16_t* data = vol.data();

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                double dx = static_cast<double>(i) - cx;
                double dy = static_cast<double>(j) - cy;
                double r2 = dx*dx + dy*dy;
                if (r2 <= outerRadius*outerRadius && r2 >= innerRadius*innerRadius) {
                    data[k * ny * nx + j * nx + i] = huValue;
                }
            }
        }
    }
}

} // anonymous namespace

// ============================================================================
// Test: Solid sphere in air → BODY contour should be generated
// ============================================================================
TEST(BodyContourGeneratorTest, SolidSphereGeneratesContours) {
    // 32x32x10 volume, sphere of HU=0 centered at (16, 16, 5) with radius 8
    auto vol = createTestVolume(32, 32, 10, -1000);
    placeSphere(*vol, 16.0, 16.0, 5.0, 8.0, 0);

    auto body = BodyContourGenerator::generate(*vol);

    ASSERT_NE(body, nullptr);
    EXPECT_EQ(body->getName(), "BODY");
    EXPECT_EQ(body->getType(), "EXTERNAL");
    EXPECT_EQ(body->getPriority(), 5);
    EXPECT_TRUE(body->isPreRasterized());
    EXPECT_GT(body->getContourCount(), 0u);
    EXPECT_GT(body->getVoxelCount(), 0u);

    // Voxel count should be much larger than contour boundary count
    // (sphere volume >> sphere surface)
    EXPECT_GT(body->getVoxelCount(), body->getContourCount() * 10);
    EXPECT_GE(body->getContourCount(), 3u);
}

// ============================================================================
// Test: All-air volume → no BODY contour generated
// ============================================================================
TEST(BodyContourGeneratorTest, AllAirReturnsNull) {
    auto vol = createTestVolume(16, 16, 4, -1000);

    auto body = BodyContourGenerator::generate(*vol);

    EXPECT_EQ(body, nullptr);
}

// ============================================================================
// Test: Entire volume is body → single contour per slice at image boundary
// ============================================================================
TEST(BodyContourGeneratorTest, FullBodyVolume) {
    auto vol = createTestVolume(16, 16, 4, 0);  // All tissue

    auto body = BodyContourGenerator::generate(*vol);

    ASSERT_NE(body, nullptr);
    // Should have one contour per slice
    EXPECT_EQ(body->getContourCount(), 4u);
}

// ============================================================================
// Test: Hollow cylinder — internal cavity should be filled
// ============================================================================
TEST(BodyContourGeneratorTest, CavityFilling) {
    // 32x32x4, hollow cylinder centered at (16,16), outer R=12, inner R=6
    auto vol = createTestVolume(32, 32, 4, -1000);
    placeHollowCylinder(*vol, 16.0, 16.0, 12.0, 6.0, 0);

    auto body = BodyContourGenerator::generate(*vol);

    ASSERT_NE(body, nullptr);
    EXPECT_GT(body->getContourCount(), 0u);

    // The contour should trace the OUTER boundary, not the inner hole.
    // Verify contour points exist and form a reasonable shape.
    const auto& contours = body->getContours();
    for (const auto& c : contours) {
        EXPECT_GE(c.points.size(), 3u) << "Each contour should have at least 3 points";
    }
}

// ============================================================================
// Test: Contour points are in patient coordinates
// ============================================================================
TEST(BodyContourGeneratorTest, ContourPointsInPatientCoords) {
    // Volume with 2mm spacing and origin at (10, 20, 30)
    auto vol = std::make_unique<Volume<int16_t>>();
    Grid grid;
    grid.setDimensions(16, 16, 4);
    grid.setSpacing(2.0, 2.0, 3.0);
    grid.setOrigin({10.0, 20.0, 30.0});
    vol->setGrid(grid);
    vol->allocate();

    // Fill with air, place a small block of tissue
    int16_t* data = vol->data();
    for (size_t i = 0; i < vol->size(); ++i) data[i] = -1000;

    // Place a 6x6 block in the center of each slice
    auto dims = grid.getDimensions();
    for (size_t k = 0; k < dims[2]; ++k) {
        for (size_t j = 5; j < 11; ++j) {
            for (size_t i = 5; i < 11; ++i) {
                data[k * dims[1] * dims[0] + j * dims[0] + i] = 0;
            }
        }
    }

    auto body = BodyContourGenerator::generate(*vol);
    ASSERT_NE(body, nullptr);

    // Check that contour z-positions reflect the patient coordinates
    // Origin Z = 30, spacing Z = 3 → slice 0 at z=30, slice 1 at z=33, etc.
    const auto& contours = body->getContours();
    for (const auto& c : contours) {
        // z should be in patient space (>= 30.0)
        EXPECT_GE(c.zPosition, 30.0);
        // All points should have the same z as zPosition
        for (const auto& pt : c.points) {
            EXPECT_DOUBLE_EQ(pt[2], c.zPosition);
        }
    }
}

// ============================================================================
// Test: Custom HU threshold
// ============================================================================
TEST(BodyContourGeneratorTest, CustomThreshold) {
    auto vol = createTestVolume(16, 16, 4, -1000);
    
    // Fill with HU = -500 (below default -400 threshold but above -600)
    int16_t* data = vol->data();
    for (size_t k = 0; k < 4; ++k) {
        for (size_t j = 4; j < 12; ++j) {
            for (size_t i = 4; i < 12; ++i) {
                data[k * 16 * 16 + j * 16 + i] = -500;
            }
        }
    }

    // Default threshold (-400) should NOT detect it
    auto bodyDefault = BodyContourGenerator::generate(*vol);
    EXPECT_EQ(bodyDefault, nullptr);

    // Custom threshold (-600) SHOULD detect it
    auto bodyCustom = BodyContourGenerator::generate(*vol, -600);
    ASSERT_NE(bodyCustom, nullptr);
    EXPECT_GT(bodyCustom->getContourCount(), 0u);
    EXPECT_GT(bodyCustom->getVoxelCount(), 0u);
}

// ============================================================================
// Test: Full-body voxel count matches expected volume
// ============================================================================
TEST(BodyContourGeneratorTest, VoxelCountMatchesVolume) {
    // 16x16x4 volume, entirely body (HU=0)
    auto vol = createTestVolume(16, 16, 4, 0);

    auto body = BodyContourGenerator::generate(*vol);
    ASSERT_NE(body, nullptr);

    // Every voxel should be in the body
    // 16 * 16 * 4 = 1024
    EXPECT_EQ(body->getVoxelCount(), 16u * 16u * 4u);
}

// ============================================================================
// Test: Pre-rasterized flag prevents re-rasterization
// ============================================================================
TEST(BodyContourGeneratorTest, SkipsReRasterization) {
    auto vol = createTestVolume(16, 16, 4, 0);

    auto body = BodyContourGenerator::generate(*vol);
    ASSERT_NE(body, nullptr);
    EXPECT_TRUE(body->isPreRasterized());

    size_t originalCount = body->getVoxelCount();
    EXPECT_GT(originalCount, 0u);

    // Calling rasterizeContours should NOT clear/recompute the indices
    Grid grid;
    grid.setDimensions(16, 16, 4);
    grid.setSpacing(1.0, 1.0, 1.0);
    grid.setOrigin({0.0, 0.0, 0.0});
    body->rasterizeContours(grid);

    EXPECT_EQ(body->getVoxelCount(), originalCount);
}

// ============================================================================
// Test: Couch/table removal — disconnected component is discarded
// ============================================================================
TEST(BodyContourGeneratorTest, CouchRemoval) {
    // 32x32x4 volume: large patient body (circle R=10 at center) + 
    // small couch strip (bottom rows 28-31, all columns) — disconnected from body
    auto vol = createTestVolume(32, 32, 4, -1000);
    auto dims = vol->getGrid().getDimensions();
    int16_t* data = vol->data();

    for (size_t k = 0; k < dims[2]; ++k) {
        size_t sliceOffset = k * dims[0] * dims[1];
        // Patient body: circle centered at (16, 16) with radius 10
        for (size_t col = 0; col < dims[1]; ++col) {
            for (size_t row = 0; row < dims[0]; ++row) {
                double dr = static_cast<double>(row) - 16.0;
                double dc = static_cast<double>(col) - 16.0;
                if (dr*dr + dc*dc <= 100.0) {  // R=10
                    data[sliceOffset + row + col * dims[0]] = 0;  // tissue
                }
            }
        }
        // Couch: flat strip at rows 28-31 across all columns (disconnected from body)
        for (size_t col = 0; col < dims[1]; ++col) {
            for (size_t row = 28; row < 32; ++row) {
                data[sliceOffset + row + col * dims[0]] = -100;  // couch material
            }
        }
    }

    auto body = BodyContourGenerator::generate(*vol);
    ASSERT_NE(body, nullptr);

    // Count the expected patient circle area (R=10): ~314 pixels per slice
    // The couch strip is 4*32 = 128 pixels per slice
    // With largest-component filtering, only the patient circle should remain
    size_t voxelsPerSlice = body->getVoxelCount() / 4;
    
    // Patient circle area ≈ π*10² ≈ 314. Should be close to that, NOT 314+128=442.
    EXPECT_GT(voxelsPerSlice, 250u);   // At least circle area (discretized)
    EXPECT_LT(voxelsPerSlice, 400u);   // Should NOT include the 128-pixel couch strip
}
