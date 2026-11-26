#include "PatientData.hpp"
#include "utils/Logger.hpp"
#include "geometry/Structure.hpp"
#include <cmath>
#include <set>

namespace optirad {

std::optional<std::array<double, 3>> PatientData::computeTargetIsocenter() const {
    if (!m_ctVolume || !m_structures) return std::nullopt;

    const auto& grid = m_ctVolume->getGrid();
    const auto dims = grid.getDimensions();
    if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) return std::nullopt;

    // Deduplicated target voxels (a voxel shared by several PTV/GTV/CTV counts once).
    std::set<size_t> targetVoxels;
    for (size_t s = 0; s < m_structures->getCount(); ++s) {
        const auto* st = m_structures->getStructure(s);
        if (!st) continue;
        const auto& type = st->getType();
        if (type == "PTV" || type == "GTV" || type == "CTV") {
            const auto& voxels = st->getVoxelIndices();
            targetVoxels.insert(voxels.begin(), voxels.end());
        }
    }
    if (targetVoxels.empty()) return std::nullopt;

    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    for (const size_t idx : targetVoxels) {
        // Column-major linear index (dims = [rows, cols, slices]): i = row, j = col, k = slice.
        const size_t k = idx / (dims[0] * dims[1]);
        const size_t j = (idx % (dims[0] * dims[1])) / dims[0];
        const size_t i = idx % dims[0];
        const Vec3 w = grid.voxelToPatient(
            {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
        sumX += w[0];
        sumY += w[1];
        sumZ += w[2];
    }

    const double n = static_cast<double>(targetVoxels.size());
    return std::array<double, 3>{sumX / n, sumY / n, sumZ / n};
}

std::optional<std::array<double, 3>> PatientData::getIsocenter() const {
    // Imported RTPLAN isocenter takes precedence and is not cached (it changes with the plan).
    if (hasImportedPlan()) {
        if (const auto iso = m_importedPlan->isocenter()) return iso;
    }

    if (!m_isocenterCache) m_isocenterCache = computeTargetIsocenter();
    return m_isocenterCache;
}

void PatientData::convertHUtoED() {
    if (!m_ctVolume) return;
    
    const auto& grid = m_ctVolume->getGrid();
    m_edVolume = std::make_unique<Volume<double>>();
    m_edVolume->setGrid(grid);
    m_edVolume->allocate();
    
    // Check that both data arrays are allocated
    if (!m_ctVolume->data() || !m_edVolume->data()) {
        Logger::error("convertHUtoED: CT or ED volume data not allocated");
        return;
    }
    
    // Simple piecewise linear HU to ED conversion
    // Based on typical HLUT values
    for (size_t i = 0; i < m_ctVolume->size(); ++i) {
        int16_t hu = m_ctVolume->data()[i];
        double ed;
        
        if (hu <= -1000) {
            ed = 0.0;  // Air
        } else if (hu <= 0) {
            // Air to water: linear interpolation
            ed = 1.0 + (hu / 1000.0);
        } else if (hu <= 100) {
            // Soft tissue
            ed = 1.0 + (hu * 0.001);
        } else {
            // Bone: steeper slope
            ed = 1.1 + (hu - 100) * 0.0005;
        }
        
        m_edVolume->data()[i] = ed;
    }
}

} // namespace optirad
