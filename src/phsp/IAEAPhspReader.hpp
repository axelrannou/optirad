#pragma once

#include "phsp/IAEAHeaderParser.hpp"
#include "phsp/PhaseSpaceData.hpp"
#include <string>
#include <vector>
#include <cstdint>

namespace optirad {

/// Reads binary IAEA phase-space files (.IAEAphsp) using the record
/// structure described in the companion .IAEAheader file.
///
/// The IAEA binary format for the TrueBeam files is:
///   - 1 sign-byte: encodes particle type via the sign of the first stored float
///                   Actually it's a combined byte where:
///                     bit pattern encodes: new_history flag (sign of charge) + particle charge
///                   For these files: photon records, electron records, positron records
///   - 6 × float32: X, Y, Z, U, V, W  (positions in cm, direction cosines unitless)
///   Total: 25 bytes per record.
///
/// The particle type is encoded in the sign of the energy value stored
/// in the binary. Since no energy is stored separately, the particle type
/// is encoded in the first byte of each record:
///   - Byte value interpretation by IAEA convention:
///     - The absolute value of the first byte encodes particle charge:
///       1 = photon, 2 = electron, 3 = positron
///     - The sign bit indicates a new history (negative = new history)
class IAEAPhspReader {
public:
    /// Read all particles from a single .IAEAphsp file.
    /// @param phspFilePath Path to the .IAEAphsp binary file
    /// @param header       Parsed header info for record structure
    /// @return PhaseSpaceData containing all particles
    static PhaseSpaceData readAll(const std::string& phspFilePath,
                                   const IAEAHeaderInfo& header);

    /// Read a subset of particles (from offset, count particles).
    /// Useful for sampling a manageable subset for visualization.
    /// @param phspFilePath Path to the .IAEAphsp binary file
    /// @param header       Parsed header info
    /// @param offset       Number of records to skip from start
    /// @param count        Maximum number of records to read
    static PhaseSpaceData readSubset(const std::string& phspFilePath,
                                      const IAEAHeaderInfo& header,
                                      int64_t offset, int64_t count);

    /// Read and concatenate particles from multiple PSF files.
    /// @param basePaths  Vector of base file paths (without extension)
    /// @param header     Parsed header info (assumed same for all files)
    /// @param maxTotal   Maximum total particles to read (0 = unlimited)
    static PhaseSpaceData readMultiple(const std::vector<std::string>& basePaths,
                                        const IAEAHeaderInfo& header,
                                        int64_t maxTotal = 0);

    /// Read an evenly-spaced sample from a single file.
    /// Reads every N-th particle to get approximately sampleSize particles.
    static PhaseSpaceData readSampled(const std::string& phspFilePath,
                                       const IAEAHeaderInfo& header,
                                       int64_t sampleSize);

private:
    /// Swap bytes of a 32-bit float for endianness conversion
    static float swapFloat(float val);

    /// Parse a single particle record from a buffer
    static Particle parseRecord(const uint8_t* buffer,
                                 const IAEAHeaderInfo& header,
                                 bool swapBytes);
};

} // namespace optirad
