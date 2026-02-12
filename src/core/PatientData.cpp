#include "PatientData.hpp"
#include "utils/Logger.hpp"
#include <cmath>

namespace optirad {

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
