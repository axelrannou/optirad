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
static constexpr uint32_t VERSION = 2; // v2: always CSR, uses direct arrays

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

    // Header
    ofs.write(MAGIC, 4);
    ofs.write(reinterpret_cast<const char*>(&VERSION), sizeof(uint32_t));

    uint64_t numVoxels   = dij.getNumVoxels();
    uint64_t numBixels   = dij.getNumBixels();
    uint64_t numNonZeros = dij.getNumNonZeros();

    ofs.write(reinterpret_cast<const char*>(&numVoxels),   sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&numBixels),   sizeof(uint64_t));
    ofs.write(reinterpret_cast<const char*>(&numNonZeros), sizeof(uint64_t));

    // Write CSR arrays directly from the matrix
    const auto& rowPtrs    = dij.getRowPtrs();
    const auto& colIndices = dij.getColIndices();
    const auto& values     = dij.getValues();

    // rowPtrs: (numVoxels+1) x uint64_t
    for (size_t i = 0; i <= numVoxels; ++i) {
        uint64_t v = rowPtrs[i];
        ofs.write(reinterpret_cast<const char*>(&v), sizeof(uint64_t));
    }
    // colIndices: nnz x uint64_t
    for (size_t i = 0; i < numNonZeros; ++i) {
        uint64_t c = colIndices[i];
        ofs.write(reinterpret_cast<const char*>(&c), sizeof(uint64_t));
    }
    // values: nnz x double
    ofs.write(reinterpret_cast<const char*>(values.data()), numNonZeros * sizeof(double));

    Logger::info("DijSerializer: Saved dij to " + filePath +
                 " (" + std::to_string(numVoxels) + "x" + std::to_string(numBixels) +
                 ", " + std::to_string(numNonZeros) + " nnz)");
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

    uint64_t numVoxels, numBixels, numNonZeros;
    ifs.read(reinterpret_cast<char*>(&numVoxels),   sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&numBixels),   sizeof(uint64_t));
    ifs.read(reinterpret_cast<char*>(&numNonZeros), sizeof(uint64_t));

    // Read CSR arrays
    std::vector<size_t> rowPtrs(numVoxels + 1);
    for (size_t i = 0; i <= numVoxels; ++i) {
        uint64_t v;
        ifs.read(reinterpret_cast<char*>(&v), sizeof(uint64_t));
        rowPtrs[i] = static_cast<size_t>(v);
    }

    std::vector<size_t> colIndices(numNonZeros);
    for (size_t i = 0; i < numNonZeros; ++i) {
        uint64_t c;
        ifs.read(reinterpret_cast<char*>(&c), sizeof(uint64_t));
        colIndices[i] = static_cast<size_t>(c);
    }

    std::vector<double> values(numNonZeros);
    ifs.read(reinterpret_cast<char*>(values.data()), numNonZeros * sizeof(double));

    // Build matrix directly from CSR (no dense allocation)
    DoseInfluenceMatrix dij;
    dij.setDimensions(numVoxels, numBixels);
    dij.loadCSR(std::move(rowPtrs), std::move(colIndices), std::move(values));

    Logger::info("DijSerializer: Loaded dij from " + filePath +
                 " (" + std::to_string(numVoxels) + "x" + std::to_string(numBixels) +
                 ", " + std::to_string(numNonZeros) + " nnz)");
    return dij;
}

bool DijSerializer::exists(const std::string& filePath) {
    return std::filesystem::exists(filePath);
}

std::string DijSerializer::buildCacheKey(
    const std::string& patientName,
    int numBeams,
    double bixelWidth,
    double doseResX)
{
    std::ostringstream oss;
    oss << patientName
        << "_" << numBeams << "beams"
        << "_bw" << std::fixed << std::setprecision(1) << bixelWidth
        << "_res" << std::fixed << std::setprecision(1) << doseResX << "mm"
        << ".dij";
    return oss.str();
}

std::string DijSerializer::getCacheDir() {
    return std::string(OPTIRAD_DATA_DIR) + "/dij_cache";
}

} // namespace optirad
