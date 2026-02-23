#include "phsp/IAEAPhspReader.hpp"
#include "utils/Logger.hpp"

#include <fstream>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace optirad {

// ---------------------------------------------------------------------------
// Byte swap helper
// ---------------------------------------------------------------------------

float IAEAPhspReader::swapFloat(float val) {
    uint32_t intVal;
    std::memcpy(&intVal, &val, 4);
    intVal = ((intVal & 0xFF000000u) >> 24) |
             ((intVal & 0x00FF0000u) >> 8)  |
             ((intVal & 0x0000FF00u) << 8)  |
             ((intVal & 0x000000FFu) << 24);
    float result;
    std::memcpy(&result, &intVal, 4);
    return result;
}

// ---------------------------------------------------------------------------
// Parse a single record
// ---------------------------------------------------------------------------

Particle IAEAPhspReader::parseRecord(const uint8_t* buffer,
                                      const IAEAHeaderInfo& header,
                                      bool swapBytes) {
    Particle p;

    // First byte: particle type encoding
    // IAEA convention: the sign byte encodes particle type
    //   |value| == 1 → photon
    //   |value| == 2 → electron
    //   |value| == 3 → positron
    // The sign bit indicates new history (negative = new primary history)
    int8_t signByte = static_cast<int8_t>(buffer[0]);
    int particleCode = std::abs(static_cast<int>(signByte));

    switch (particleCode) {
        case 1: p.type = ParticleType::Photon;   break;
        case 2: p.type = ParticleType::Electron;  break;
        case 3: p.type = ParticleType::Positron;  break;
        default: p.type = ParticleType::Photon;   break; // fallback
    }

    // Weight is constant for these files
    p.weight = header.constantWeight;

    // Read float fields starting at offset 1
    size_t offset = 1;
    auto readFloat = [&]() -> float {
        float val;
        std::memcpy(&val, buffer + offset, 4);
        offset += 4;
        if (swapBytes) val = swapFloat(val);
        return val;
    };

    // IAEA binary field order: sign, [Ek], [X], [Y], [Z], [U], [V], [Wt], [extras]
    // Energy sign encodes W direction: negative energy → W is negative.
    float rawEnergy = 0.0f;
    if (header.storesEnergy) rawEnergy = readFloat();

    if (header.storesX) p.position[0] = static_cast<double>(readFloat()) * 10.0; // cm → mm
    if (header.storesY) p.position[1] = static_cast<double>(readFloat()) * 10.0;
    if (header.storesZ) p.position[2] = static_cast<double>(readFloat()) * 10.0; // scoring plane Z

    float u = 0.0f, v = 0.0f;
    if (header.storesU) u = readFloat();
    if (header.storesV) v = readFloat();
    p.direction[0] = static_cast<double>(u);
    p.direction[1] = static_cast<double>(v);

    // Compute W = sqrt(1 - U² - V²), sign from energy sign
    double uv2 = static_cast<double>(u) * u + static_cast<double>(v) * v;
    double w = (uv2 < 1.0) ? std::sqrt(1.0 - uv2) : 0.0;
    if (rawEnergy < 0.0f) w = -w;
    p.direction[2] = w;

    if (header.storesWeight) {
        p.weight = static_cast<double>(readFloat());
    }

    // Energy: use absolute value (sign was for W direction)
    p.energy = static_cast<double>(std::fabs(rawEnergy));

    // If energy not stored, fall back to mean from header
    if (!header.storesEnergy) {
        switch (p.type) {
            case ParticleType::Photon:   p.energy = header.photonMeanEnergy;   break;
            case ParticleType::Electron: p.energy = header.electronMeanEnergy; break;
            case ParticleType::Positron: p.energy = header.positronMeanEnergy; break;
        }
    }

    return p;
}

// ---------------------------------------------------------------------------
// Read all particles
// ---------------------------------------------------------------------------

PhaseSpaceData IAEAPhspReader::readAll(const std::string& phspFilePath,
                                        const IAEAHeaderInfo& header) {
    return readSubset(phspFilePath, header, 0, header.totalParticles);
}

// ---------------------------------------------------------------------------
// Read subset
// ---------------------------------------------------------------------------

PhaseSpaceData IAEAPhspReader::readSubset(const std::string& phspFilePath,
                                            const IAEAHeaderInfo& header,
                                            int64_t offset, int64_t count) {
    Logger::debug("Reading IAEA phsp: " + phspFilePath +
                  " (offset=" + std::to_string(offset) +
                  ", count=" + std::to_string(count) + ")");

    std::ifstream file(phspFilePath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("IAEAPhspReader: cannot open file: " + phspFilePath);
    }

    bool swapBytes = header.needsByteSwap();
    int recLen = header.recordLength;

    // Seek to offset
    if (offset > 0) {
        file.seekg(static_cast<std::streamoff>(offset) * recLen);
        if (!file.good()) {
            throw std::runtime_error("IAEAPhspReader: seek failed in " + phspFilePath);
        }
    }

    PhaseSpaceData data;
    data.reserve(static_cast<size_t>(std::min(count, int64_t(10000000)))); // cap at 10M for reservation

    std::vector<uint8_t> buffer(recLen);
    int64_t read = 0;

    while (read < count && file.read(reinterpret_cast<char*>(buffer.data()), recLen)) {
        Particle p = parseRecord(buffer.data(), header, swapBytes);
        data.addParticle(p);
        ++read;
    }

    Logger::debug("Read " + std::to_string(read) + " particles from " + phspFilePath);
    return data;
}

// ---------------------------------------------------------------------------
// Read sampled (every N-th record)
// ---------------------------------------------------------------------------

PhaseSpaceData IAEAPhspReader::readSampled(const std::string& phspFilePath,
                                             const IAEAHeaderInfo& header,
                                             int64_t sampleSize) {
    Logger::debug("Sampling ~" + std::to_string(sampleSize) + " particles from " + phspFilePath);

    std::ifstream file(phspFilePath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("IAEAPhspReader: cannot open file: " + phspFilePath);
    }

    bool swapBytes = header.needsByteSwap();
    int recLen = header.recordLength;

    // Compute stride: read every Nth record
    int64_t totalRecords = header.totalParticles;
    int64_t stride = std::max(int64_t(1), totalRecords / sampleSize);

    PhaseSpaceData data;
    data.reserve(static_cast<size_t>(sampleSize + 100));

    std::vector<uint8_t> buffer(recLen);
    int64_t recordIndex = 0;
    int64_t read = 0;

    while (read < sampleSize && file.read(reinterpret_cast<char*>(buffer.data()), recLen)) {
        if (recordIndex % stride == 0) {
            Particle p = parseRecord(buffer.data(), header, swapBytes);
            data.addParticle(p);
            ++read;
        }
        ++recordIndex;

        // If stride > 1, skip ahead by seeking
        if (stride > 1 && read < sampleSize) {
            int64_t skipRecords = stride - 1;
            file.seekg(static_cast<std::streamoff>(skipRecords) * recLen, std::ios::cur);
            recordIndex += skipRecords;
        }
    }

    Logger::debug("Sampled " + std::to_string(read) + " particles (stride=" +
                  std::to_string(stride) + ")");
    return data;
}

// ---------------------------------------------------------------------------
// Read multiple files
// ---------------------------------------------------------------------------

PhaseSpaceData IAEAPhspReader::readMultiple(const std::vector<std::string>& basePaths,
                                              const IAEAHeaderInfo& header,
                                              int64_t maxTotal) {
    PhaseSpaceData combined;
    int64_t totalRead = 0;

    for (const auto& basePath : basePaths) {
        std::string phspPath = basePath + ".IAEAphsp";

        int64_t remaining = (maxTotal > 0) ? (maxTotal - totalRead) : header.totalParticles;
        if (remaining <= 0) break;

        PhaseSpaceData chunk = readSubset(phspPath, header, 0, remaining);
        totalRead += static_cast<int64_t>(chunk.size());

        // Merge into combined
        auto& dst = combined.particles();
        auto& src = chunk.particles();
        dst.insert(dst.end(),
                   std::make_move_iterator(src.begin()),
                   std::make_move_iterator(src.end()));

        if (maxTotal > 0 && totalRead >= maxTotal) break;
    }

    Logger::debug("Read " + std::to_string(totalRead) + " particles from " +
                  std::to_string(basePaths.size()) + " files");
    return combined;
}

} // namespace optirad
