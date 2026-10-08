#pragma once

#include "DoseInfluenceMatrix.hpp"
#include <string>

namespace optirad {

/**
 * Serializes/deserializes DoseInfluenceMatrix to/from binary files.
 * 
 * Binary format (v4):
 *   - Magic bytes: "ODIJ" (4 bytes)
 *   - Version: uint32_t (3)
 *   - numVoxels, numBixels, numBlocks: uint64_t
 *   - per block: nnz (uint64_t), rowPtrs (numVoxels+1) x uint32_t,
 *                colIndices nnz x uint32_t, values nnz x float
 *   - computed-row mask: size (uint64_t, 0 = all rows) followed by size x uint8_t
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
     * @param relativeThreshold  Dij relative threshold (fraction)
     * @param excludeExternal    Body-only voxels were skipped
     * @return Filename like "JOHN_DOE_90beams_bw5.0_res2.5mm_thr1e-04_noext_e5.dij"
     */
    static std::string buildCacheKey(
        const std::string& patientName,
        int numBeams,
        double bixelWidth,
        double doseResX,
        double relativeThreshold = 0.0,
        bool excludeExternal = false);

    /// Bump whenever the dose engine changes the Dij values, so stale caches are not reused.
    static constexpr int kEngineVersion = 5;

    /**
     * Get the default cache directory path.
     */
    static std::string getCacheDir();
};

} // namespace optirad
