#include <gtest/gtest.h>
#include "dose/DijSerializer.hpp"
#include "dose/DoseInfluenceMatrix.hpp"

#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <cstring>

namespace optirad::tests {

class DijSerializerTest : public ::testing::Test {
protected:
    std::string m_tmpDir;

    void SetUp() override {
        m_tmpDir = std::filesystem::temp_directory_path() / "optirad_test_dij";
        std::filesystem::create_directories(m_tmpDir);
    }

    void TearDown() override {
        std::filesystem::remove_all(m_tmpDir);
    }

    std::string tmpPath(const std::string& name) const {
        return m_tmpDir + "/" + name;
    }
};

// ============================================================================
// Round-trip save → load
// ============================================================================

TEST_F(DijSerializerTest, SaveLoadRoundTrip) {
    DoseInfluenceMatrix dij(10, 5);
    dij.setValue(0, 0, 1.5);
    dij.setValue(0, 3, 2.5);
    dij.setValue(3, 1, 3.5);
    dij.setValue(9, 4, 4.5);
    dij.finalize();

    std::string path = tmpPath("test.dij");
    ASSERT_TRUE(DijSerializer::save(dij, path));
    EXPECT_TRUE(DijSerializer::exists(path));

    auto loaded = DijSerializer::load(path);
    EXPECT_EQ(loaded.getNumVoxels(), 10u);
    EXPECT_EQ(loaded.getNumBixels(), 5u);
    EXPECT_EQ(loaded.getNumNonZeros(), 4u);
    EXPECT_TRUE(loaded.isFinalized());

    EXPECT_DOUBLE_EQ(loaded.getValue(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(loaded.getValue(0, 3), 2.5);
    EXPECT_DOUBLE_EQ(loaded.getValue(3, 1), 3.5);
    EXPECT_DOUBLE_EQ(loaded.getValue(9, 4), 4.5);
    EXPECT_DOUBLE_EQ(loaded.getValue(5, 2), 0.0);
}

TEST_F(DijSerializerTest, SaveLoadComputeDoseConsistent) {
    DoseInfluenceMatrix dij(5, 3);
    dij.setValue(0, 0, 1.0);
    dij.setValue(0, 1, 2.0);
    dij.setValue(1, 2, 3.0);
    dij.setValue(4, 0, 4.0);
    dij.finalize();

    std::vector<double> weights = {1.0, 2.0, 3.0};
    auto doseBefore = dij.computeDose(weights);

    std::string path = tmpPath("dose_check.dij");
    ASSERT_TRUE(DijSerializer::save(dij, path));
    auto loaded = DijSerializer::load(path);
    auto doseAfter = loaded.computeDose(weights);

    ASSERT_EQ(doseBefore.size(), doseAfter.size());
    for (size_t i = 0; i < doseBefore.size(); ++i) {
        EXPECT_DOUBLE_EQ(doseBefore[i], doseAfter[i]) << "mismatch at voxel " << i;
    }
}

TEST_F(DijSerializerTest, EmptyMatrixRoundTrip) {
    DoseInfluenceMatrix dij(100, 50);
    dij.finalize();

    std::string path = tmpPath("empty.dij");
    ASSERT_TRUE(DijSerializer::save(dij, path));

    auto loaded = DijSerializer::load(path);
    EXPECT_EQ(loaded.getNumVoxels(), 100u);
    EXPECT_EQ(loaded.getNumBixels(), 50u);
    EXPECT_EQ(loaded.getNumNonZeros(), 0u);
}

// ============================================================================
// Invalid file handling
// ============================================================================

TEST_F(DijSerializerTest, LoadNonExistentThrows) {
    EXPECT_THROW(DijSerializer::load(tmpPath("nonexistent.dij")), std::runtime_error);
}

TEST_F(DijSerializerTest, LoadCorruptMagicThrows) {
    std::string path = tmpPath("corrupt.dij");
    {
        std::ofstream ofs(path, std::ios::binary);
        ofs.write("XXXX", 4); // wrong magic
    }
    EXPECT_THROW(DijSerializer::load(path), std::runtime_error);
}

// ============================================================================
// exists()
// ============================================================================

TEST_F(DijSerializerTest, ExistsReturnsFalseForMissing) {
    EXPECT_FALSE(DijSerializer::exists(tmpPath("nope.dij")));
}

// ============================================================================
// buildCacheKey
// ============================================================================

TEST_F(DijSerializerTest, BuildCacheKeyDeterministic) {
    auto key1 = DijSerializer::buildCacheKey("TEST", 5, 5.0, 2.5);
    auto key2 = DijSerializer::buildCacheKey("TEST", 5, 5.0, 2.5);
    EXPECT_EQ(key1, key2);

    // Contains expected components
    EXPECT_NE(key1.find("TEST"), std::string::npos);
    EXPECT_NE(key1.find("5beams"), std::string::npos);
    EXPECT_NE(key1.find(".dij"), std::string::npos);
}

TEST_F(DijSerializerTest, BuildCacheKeyVariesWithParams) {
    auto k1 = DijSerializer::buildCacheKey("A", 1, 5.0, 2.5);
    auto k2 = DijSerializer::buildCacheKey("B", 1, 5.0, 2.5);
    auto k3 = DijSerializer::buildCacheKey("A", 2, 5.0, 2.5);
    auto k4 = DijSerializer::buildCacheKey("A", 1, 10.0, 2.5);
    EXPECT_NE(k1, k2);
    EXPECT_NE(k1, k3);
    EXPECT_NE(k1, k4);
}

} // namespace optirad::tests
