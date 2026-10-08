#include "DijSerializer.hpp"
#include "utils/Logger.hpp"

#include <fstream>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <iomanip>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad {

static constexpr char MAGIC[4] = {'O', 'D', 'I', 'J'};
static constexpr uint32_t VERSION = 4; // v4: adds the computed-row mask to v3 (CSR column blocks)

namespace {

template <typename T>
void writeVec(std::ofstream& ofs, const std::vector<T>& v) {
    ofs.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
}

template <typename T>
void readVec(std::ifstream& ifs, std::vector<T>& v, size_t n) {
    v.resize(n);
    ifs.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
}

} // namespace

bool DijSerializer::save(const DoseInfluenceMatrix& dij, const std::string& filePath) {
    if (!dij.isFinalized()) {
        Logger::error("DijSerializer: Matrix must be finalized before saving.");
        return false;
    }

    auto parentDir = std::filesystem::path(filePath).parent_path();
    if (!parentDir.empty()) {
        std::filesystem::create_directories(parentDir);
    }

    std::ofstream ofs(filePath, std::ios::binary);
    if (!ofs) {
        Logger::error("DijSerializer: Cannot open file for writing: " + filePath);
        return false;
    }

    // Header: magic, version, numVoxels, numBixels, numBlocks
    ofs.write(MAGIC, 4);
    ofs.write(reinterpret_cast<const char*>(&VERSION), sizeof(uint32_t));

    const uint64_t numVoxels = dij.getNumVoxels();
    const uint64_t numBixels = dij.getNumBixels();
    const uint64_t numBlocks = dij.getBlocks().size();
    ofs.write(reinterpret_cast<const char*>(&numVoxels), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&numBixels), sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&numBlocks), sizeof(uint64_t));

    for (const auto& block : dij.getBlocks()) {
        const uint64_t nnz = block.values.size();
        ofs.write(reinterpret_cast<const char*>(&nnz), sizeof(uint64_t));
        writeVec(ofs, block.rowPtrs);
        writeVec(ofs, block.colIndices);
        writeVec(ofs, block.values);
    }

    const uint64_t maskSize = dij.getRowMask().size();
    ofs.write(reinterpret_cast<const char*>(&maskSize), sizeof(uint64_t));
    writeVec(ofs, dij.getRowMask());

    if (!ofs) {
        Logger::error("DijSerializer: Write failed (disk full?): " + filePath);
        return false;
    }

    Logger::info("DijSerializer: Saved dij to " + filePath +
                 " (" + std::to_string(numVoxels) + "x" + std::to_string(numBixels) +
                 ", " + std::to_string(dij.getNumNonZeros()) + " nnz, " +
                 std::to_string(numBlocks) + " blocks)");
    return true;
}

DoseInfluenceMatrix DijSerializer::load(const std::string& filePath) {
    std::ifstream ifs(filePath, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("DijSerializer: Cannot open file: " + filePath);
    }

    char magic[4];
    ifs.read(magic, 4);
    if (magic[0] != MAGIC[0] || magic[1] != MAGIC[1] ||
        magic[2] != MAGIC[2] || magic[3] != MAGIC[3]) {
        throw std::runtime_error("DijSerializer: Invalid magic bytes in " + filePath);
    }

    uint32_t version;
    ifs.read(reinterpret_cast<char*>(&version), sizeof(uint32_t));
    if (version != VERSION) {
        throw std::runtime_error("DijSerializer: Unsupported version " +
                                 std::to_string(version) + " (expected " +
                                 std::to_string(VERSION) + ")");
    }

    uint64_t numVoxels, numBixels, numBlocks;
    ifs.read(reinterpret_cast<char*>(&numVoxels), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&numBixels), sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&numBlocks), sizeof(uint64_t));

    DoseInfluenceMatrix dij;
    dij.setDimensions(numVoxels, numBixels);

    for (uint64_t b = 0; b < numBlocks; ++b) {
        uint64_t nnz;
        ifs.read(reinterpret_cast<char*>(&nnz), sizeof(uint64_t));
        DoseInfluenceMatrix::Block block;
        readVec(ifs, block.rowPtrs, numVoxels + 1);
        readVec(ifs, block.colIndices, nnz);
        readVec(ifs, block.values, nnz);
        if (!ifs) throw std::runtime_error("DijSerializer: Truncated file: " + filePath);
        dij.addBlock(std::move(block));
    }

    uint64_t maskSize = 0;
    ifs.read(reinterpret_cast<char*>(&maskSize), sizeof(uint64_t));
    if (maskSize > 0) {
        std::vector<uint8_t> mask;
        readVec(ifs, mask, maskSize);
        if (!ifs) throw std::runtime_error("DijSerializer: Truncated row mask: " + filePath);
        dij.setComputedRows(std::move(mask));
    }

    Logger::info("DijSerializer: Loaded dij from " + filePath +
                 " (" + std::to_string(numVoxels) + "x" + std::to_string(numBixels) +
                 ", " + std::to_string(dij.getNumNonZeros()) + " nnz)");
    return dij;
}

bool DijSerializer::exists(const std::string& filePath) {
    return std::filesystem::exists(filePath);
}

std::string DijSerializer::buildCacheKey(
    const std::string& patientName,
    int numBeams,
    double bixelWidth,
    double doseResX,
    double relativeThreshold,
    bool excludeExternal)
{
    std::ostringstream oss;
    oss << patientName
        << "_" << numBeams << "beams"
        << "_bw" << std::fixed << std::setprecision(1) << bixelWidth
        << "_res" << std::fixed << std::setprecision(1) << doseResX << "mm"
        << "_thr" << std::scientific << std::setprecision(0) << relativeThreshold
        << (excludeExternal ? "_noext" : "")
        << "_e" << kEngineVersion
        << ".dij";
    return oss.str();
}

std::string DijSerializer::getCacheDir() {
    return std::string(OPTIRAD_DATA_DIR) + "/dij_cache";
}

} // namespace optirad
