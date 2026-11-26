#include "Plan.hpp"
#include "PatientData.hpp"
#include "geometry/Volume.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <numeric>
#include <vector>
#include <set>

namespace optirad {

void Plan::setName(const std::string& name) { m_name = name; }

double Plan::machineMuPerWeight() const {
    const double sad = m_machine.getSAD() > 0.0 ? m_machine.getSAD() : 1000.0;
    const int fractions = m_numOfFractions > 0 ? m_numOfFractions : 1;
    return 100.0 * (sad / 1000.0) * (sad / 1000.0) / fractions;
}
const std::string& Plan::getName() const { return m_name; }

void Plan::addBeam(const Beam& beam) { m_beams.push_back(beam); }
const std::vector<Beam>& Plan::getBeams() const { return m_beams; }
size_t Plan::getNumBeams() const { return m_beams.size(); }

std::array<double, 3> Plan::computeIsoCenter() const {
    if (!m_patientData) {
        return {0.0, 0.0, 0.0};
    }

    auto* ct = m_patientData->getCTVolume();
    auto* structureSet = m_patientData->getStructureSet();
    
    if (!ct || !structureSet) {
        return {0.0, 0.0, 0.0};
    }

    // Report diagnostics before delegating to the shared computation in PatientData.
    int ptvCount = 0;
    for (size_t i = 0; i < structureSet->getCount(); ++i) {
        const auto* structure = structureSet->getStructure(i);
        // I will fix it later but to compare with matRad we keep all "targets" type
        if (structure && structure->getType() == "PTV" || structure->getType() == "GTV" || structure->getType() == "CTV") {
            ptvCount++;
            if (structure->getVoxelIndices().empty()) {
                Logger::warn("PTV structure '" + structure->getName() + 
                           "' has no voxel indices (contours not rasterized yet)");
            }
        }
    }

    // Error out if no PTV structures found
    if (ptvCount == 0) {
        Logger::error("No PTV structures found in structure set");
        Logger::error("Available structures:");
        for (size_t i = 0; i < structureSet->getCount(); ++i) {
            const auto* s = structureSet->getStructure(i);
            if (s) {
                Logger::error("  - " + s->getName() + " (type: " + s->getType() + ")");
            }
        }
        Logger::error("Cannot compute isocenter without PTV structures");
        return {0.0, 0.0, 0.0};
    }

    const auto iso = m_patientData->computeTargetIsocenter();
    if (!iso) {
        Logger::error("PTV structures found but have no voxel indices");
        Logger::error("You must rasterize structure contours before computing isocenter");
        return {0.0, 0.0, 0.0};
    }

    Logger::info("Computed isocenter from target structures");

    return *iso;
}

void Plan::printSummary() const {
    std::cout << "\n=== Treatment Plan Summary ===\n";
    std::cout << "  Plan Name:       " << (m_name.empty() ? "(unnamed)" : m_name) << "\n";
    std::cout << "  Radiation Mode:  " << m_radiationMode << "\n";
    std::cout << "  Machine:         " << m_machine.getName() << "\n";
    std::cout << "  SAD:             " << m_machine.getSAD() << " mm\n";
    std::cout << "  SCD:             " << m_machine.getSCD() << " mm\n";
    std::cout << "  Energy:          " << m_machine.getData().energy << " MV\n";
    std::cout << "  Num Fractions:   " << m_numOfFractions << "\n";
    if (m_deliveryMode == DeliveryMode::VMAT) {
        std::cout << "  Delivery:        VMAT (" << m_vmat.arcs.size() << " arc(s), "
                  << m_vmat.optimizationSpacingDeg << " deg optimization beams, "
                  << m_vmat.controlPointSpacingDeg << " deg control points)\n";
        for (size_t a = 0; a < m_vmat.arcs.size(); ++a) {
            const auto& arc = m_vmat.arcs[a];
            std::cout << "    Arc " << a + 1 << ": " << arc.startAngle << " -> " << arc.stopAngle
                      << (arc.clockwise ? " CW" : " CCW") << ", couch " << arc.couchAngle
                      << ", collimator " << arc.collimatorAngle << "\n";
        }
    } else {
        std::cout << "  Delivery:        SMLC (step-and-shoot)\n";
    }

    std::cout << "\n  STF Properties:\n";
    std::cout << "    Num Beams:     " << m_stfProperties.numOfBeams << "\n";
    std::cout << "    Bixel Width:   " << m_stfProperties.bixelWidth << " mm\n";

    if (!m_stfProperties.gantryAngles.empty()) {
        std::cout << "    Gantry Angles: [" << m_stfProperties.gantryAngles.front()
                  << " : ... : " << m_stfProperties.gantryAngles.back() << "] deg"
                  << " (" << m_stfProperties.gantryAngles.size() << " angles)\n";

        if (m_stfProperties.gantryAngles.size() > 1) {
            double step = m_stfProperties.gantryAngles[1] - m_stfProperties.gantryAngles[0];
            std::cout << "    Gantry Step:   " << step << " deg\n";
        }
    }

    if (!m_stfProperties.isoCenters.empty()) {
        const auto& iso = m_stfProperties.isoCenters[0];
        std::cout << "    Isocenter:     [" << iso[0] << ", " << iso[1] << ", " << iso[2] << "] mm\n";
    }

    // Machine constraints
    const auto& c = m_machine.getConstraints();
    std::cout << "\n  Machine Constraints:\n";
    std::cout << "    Gantry Speed:  [" << c.gantryRotationSpeed[0] << ", " << c.gantryRotationSpeed[1] << "] deg/s\n";
    std::cout << "    Leaf Speed:    [" << c.leafSpeed[0] << ", " << c.leafSpeed[1] << "] mm/s\n";
    std::cout << "    MU Rate:       [" << c.monitorUnitRate[0] << ", " << c.monitorUnitRate[1] << "] MU/s\n";

    std::cout << "==============================\n";
}

} // namespace optirad
