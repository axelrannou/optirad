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
    
    if (!m_beams.empty()) {
        std::cout << "\nBeam Details:\n";
        for (size_t i = 0; i < m_beams.size(); ++i) {
            const auto& beam = m_beams[i];
            std::cout << "  Beam " << i << ":\n";
            std::cout << "    Gantry:      " << beam.getGantryAngle() << " deg\n";
            std::cout << "    Couch:       " << beam.getCouchAngle() << " deg\n";
            std::cout << "    Isocenter:   [" << beam.getIsocenter()[0] << ", " 
                      << beam.getIsocenter()[1] << ", " << beam.getIsocenter()[2] << "] mm\n";
            std::cout << "    Rays:        " << beam.getNumOfRays() << "\n";
            std::cout << "    Bixels:      " << beam.getTotalNumOfBixels() << "\n";
            std::cout << "    Bixel Width: " << beam.getBixelWidth() << " mm\n";
            std::cout << "    SAD:         " << beam.getSAD() << " mm\n";
        }
    }
    std::cout << "===================\n";
}

} // namespace optirad
