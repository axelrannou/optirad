#pragma once

#include <string>
#include <array>
#include <cstdint>

namespace optirad {

/// Information parsed from an IAEA phase-space header file (.IAEAheader).
/// See IAEA(NDS)-0484 technical report for the specification.
struct IAEAHeaderInfo {
    // File identification
    int iaeaIndex = 0;
    std::string title;
    int fileType = 0; // 0 = event generator

    // Record structure — IAEA RECORD_CONTENTS indices:
    //   [0]=Ek, [1]=X, [2]=Y, [3]=Z, [4]=U, [5]=V, [6]=Wt, [7]=ExF, [8]=ExL
    // Note: W (direction cosine Z) is NEVER stored; it is computed from U,V.
    // The sign of W is encoded in the sign of the Energy float.
    bool storesEnergy = false;
    bool storesX = false, storesY = false, storesZ = false;
    bool storesU = false, storesV = false;
    bool storesWeight = false;
    int numExtraFloats = 0;
    int numExtraLongs = 0;
    double constantWeight = 1.0;
    int recordLength = 25; // bytes per particle record
    int byteOrder = 1234;  // 1234 = little-endian, 4321 = big-endian

    // Particle counts
    int64_t origHistories = 0;
    int64_t totalParticles = 0;
    int64_t numPhotons = 0;
    int64_t numElectrons = 0;
    int64_t numPositrons = 0;

    // Energy info
    double photonEnergyCutoff = 0.0;   // MeV
    double particleEnergyCutoff = 0.0; // MeV

    // Geometry extents (cm)
    std::array<double, 2> xRange = {0.0, 0.0};
    std::array<double, 2> yRange = {0.0, 0.0};
    std::array<double, 2> zRange = {0.0, 0.0};

    // Statistical info (mean energy in MeV)
    double photonMeanEnergy = 0.0;
    double electronMeanEnergy = 0.0;
    double positronMeanEnergy = 0.0;

    // MC code
    std::string mcCodeVersion;

    // Nominal SSD (cm)
    double nominalSSD = 100.0;

    // Coordinate system description
    std::string coordinateSystemDesc;

    /// Compute expected record length from stored fields.
    int computeRecordLength() const {
        // 1 byte for particle type sign + 4 bytes per stored float field
        int nFloats = (storesEnergy ? 1 : 0)
                    + (storesX ? 1 : 0) + (storesY ? 1 : 0) + (storesZ ? 1 : 0)
                    + (storesU ? 1 : 0) + (storesV ? 1 : 0)
                    + (storesWeight ? 1 : 0) + numExtraFloats;
        return 1 + nFloats * 4 + numExtraLongs * 4;
    }

    /// Check if host byte order matches file byte order
    bool needsByteSwap() const;
};

/// Parses IAEA phase-space header files (.IAEAheader).
class IAEAHeaderParser {
public:
    /// Parse a header file and return the extracted information.
    /// @throws std::runtime_error on file open or parse errors
    static IAEAHeaderInfo parse(const std::string& headerFilePath);
};

} // namespace optirad
