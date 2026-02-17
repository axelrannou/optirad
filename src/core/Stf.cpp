#include "Stf.hpp"
#include "../utils/Logger.hpp"
#include <iostream>

namespace optirad {

void Stf::addBeam(const Beam& beam) {
    m_beams.push_back(beam);
}

void Stf::addBeam(Beam&& beam) {
    m_beams.push_back(std::move(beam));
}

const Beam* Stf::getBeam(size_t index) const {
    if (index >= m_beams.size()) {
        return nullptr;
    }
    return &m_beams[index];
}

Beam* Stf::getBeam(size_t index) {
    if (index >= m_beams.size()) {
        return nullptr;
    }
    return &m_beams[index];
}

size_t Stf::getTotalNumOfRays() const {
    size_t total = 0;
    for (const auto& beam : m_beams) {
        total += beam.getNumOfRays();
    }
    return total;
}

size_t Stf::getTotalNumOfBixels() const {
    size_t total = 0;
    for (const auto& beam : m_beams) {
        total += beam.getTotalNumOfBixels();
    }
    return total;
}

void Stf::computeAllSourcePoints() {
    for (auto& beam : m_beams) {
        beam.computeSourcePoints();
    }
}

void Stf::generateAllRays(double bixelWidth, const std::array<double, 2>& fieldSize) {
    for (auto& beam : m_beams) {
        beam.generateRays(bixelWidth, fieldSize);
    }
}

void Stf::printSummary() const {
    std::cout << "\n=== STF Summary ===\n";
    std::cout << "Number of beams: " << getCount() << "\n";
    std::cout << "Total rays:      " << getTotalNumOfRays() << "\n";
    std::cout << "Total bixels:    " << getTotalNumOfBixels() << "\n";

    #ifndef NDEBUG
    // Show first 10 rays in Debug builds
    if (!m_beams.empty()) {
        std::cout << "\nBeam Details:\n";
        // only first 5 beams for brevity
        for (size_t i = 0; i < m_beams.size() && i < 5; ++i) {
            const auto& beam = m_beams[i];
            std::cout << "  Beam " << i << ":\n";
            std::cout << "    Gantry:        " << beam.getGantryAngle() << " deg\n";
            std::cout << "    Couch:         " << beam.getCouchAngle() << " deg\n";
            std::cout << "    Isocenter:     [" << beam.getIsocenter()[0] << ", " 
                      << beam.getIsocenter()[1] << ", " << beam.getIsocenter()[2] << "] mm\n";
            std::cout << "    Bixel Width:   " << beam.getBixelWidth() << " mm\n";
            std::cout << "    Radiation:     " << beam.getRadiationMode() << "\n";
            std::cout << "    Machine:       " << beam.getMachineName() << "\n";
            std::cout << "    SAD:           " << beam.getSAD() << " mm\n";
            std::cout << "    SCD:           " << beam.getSCD() << " mm\n";
            std::cout << "    SourcePt BEV:  [" << beam.getSourcePointBev()[0] << ", "
                      << beam.getSourcePointBev()[1] << ", " << beam.getSourcePointBev()[2] << "]\n";
            std::cout << "    SourcePt LPS:  [" << beam.getSourcePoint()[0] << ", "
                      << beam.getSourcePoint()[1] << ", " << beam.getSourcePoint()[2] << "]\n";
            std::cout << "    Rays:          " << beam.getNumOfRays() << "\n";
            std::cout << "    Total bixels:  " << beam.getTotalNumOfBixels() << "\n";
            
            size_t numRaysToPrint = std::min(static_cast<size_t>(10), beam.getNumOfRays());
            for (size_t r = 0; r < numRaysToPrint; ++r) {
                const auto* ray = beam.getRay(r);
                std::cout << "    Ray " << r << ":\n";
                std::cout << "      rayPos_bev:    [" << ray->getRayPosBev()[0] << ", "
                          << ray->getRayPosBev()[1] << ", " << ray->getRayPosBev()[2] << "]\n";
                std::cout << "      targetPt_bev:  [" << ray->getTargetPointBev()[0] << ", "
                          << ray->getTargetPointBev()[1] << ", " << ray->getTargetPointBev()[2] << "]\n";
                std::cout << "      rayPos:        [" << ray->getRayPos()[0] << ", "
                          << ray->getRayPos()[1] << ", " << ray->getRayPos()[2] << "]\n";
                std::cout << "      targetPt:      [" << ray->getTargetPoint()[0] << ", "
                          << ray->getTargetPoint()[1] << ", " << ray->getTargetPoint()[2] << "]\n";
                std::cout << "      energy:        " << ray->getEnergy() << "\n";
                
                const auto& cornersIso = ray->getBeamletCornersAtIso();
                std::cout << "      beamletCornersAtIso:\n";
                for (size_t c = 0; c < 4; ++c) {
                    std::cout << "        [" << cornersIso[c][0] << ", " 
                              << cornersIso[c][1] << ", " << cornersIso[c][2] << "]\n";
                }
                
                const auto& cornersSCD = ray->getRayCornersSCD();
                std::cout << "      rayCornersSCD:\n";
                for (size_t c = 0; c < 4; ++c) {
                    std::cout << "        [" << cornersSCD[c][0] << ", " 
                              << cornersSCD[c][1] << ", " << cornersSCD[c][2] << "]\n";
                }
            }
        }
    }
    #endif
    std::cout << "===================\n";
}

} // namespace optirad
