#pragma once

#include "DoseInfluenceMatrix.hpp"
#include <string>

namespace optirad {

/**
 * Serializes/deserializes DoseInfluenceMatrix to/from binary files.
 * 
 * Binary format (v2):
 *   - Magic bytes: "ODIJ" (4 bytes)
 *   - Version: uint32_t (2)
 *   - numVoxels: uint64_t
 *   - numBixels: uint64_t
 *   - numNonZeros: uint64_t
 *   - rowPtrs: (numVoxels+1) x uint64_t
 *   - colIndices: numNonZeros x uint64_t
 *   - values: numNonZeros x double
 */
class DijSerializer {
public:
    /**
     * Save a DoseInfluenceMatrix to a binary file.
     * The matrix should be finalized before saving.
     */
    static bool save(const DoseInfluenceMatrix& dij, const std::string& filePath);

    /**
     * Load a DoseInfluenceMatrix from a binary file.
     */
    static DoseInfluenceMatrix load(const std::string& filePath);

    /**
     * Check if a cache file exists.
     */
    static bool exists(const std::string& filePath);

    /**
     * Build a deterministic cache filename.
     * @param patientName   Patient identifier
     * @param numBeams      Number of beams
     * @param bixelWidth    Bixel width in mm
     * @param doseResX      Dose grid resolution x (mm)
     * @return Filename like "JOHN_DOE_90beams_bw5.0_res2.5mm.dij"
     */
    static std::string buildCacheKey(
        const std::string& patientName,
        int numBeams,
        double bixelWidth,
        double doseResX);

    /**
     * Get the default cache directory path.
     */
    static std::string getCacheDir();
};

} // namespace optirad
