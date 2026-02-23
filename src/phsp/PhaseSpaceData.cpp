#include "phsp/PhaseSpaceData.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace optirad {

PhaseSpaceData PhaseSpaceData::sample(size_t n, unsigned int seed) const {
    PhaseSpaceData result;
    if (m_particles.empty() || n == 0) return result;

    if (n >= m_particles.size()) {
        result.m_particles = m_particles;
        return result;
    }

    // Create index array and shuffle
    std::vector<size_t> indices(m_particles.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);

    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        result.addParticle(m_particles[indices[i]]);
    }

    return result;
}

void PhaseSpaceData::filterByJaws(double jawX1, double jawX2,
                                   double jawY1, double jawY2,
                                   double SAD, double scoringPlaneZ) {
    // Project each particle position to the isocenter plane
    // The scoring plane is at some Z position above the isocenter; 
    // we project along the particle direction to Z = 0 (isocenter plane)
    m_particles.erase(
        std::remove_if(m_particles.begin(), m_particles.end(),
            [&](const Particle& p) {
                // Direction cosines: u, v, w (w is along beam axis toward patient)
                double w = p.direction[2];
                if (std::abs(w) < 1e-10) return true; // particle going sideways, remove

                // Project from scoring plane position to isocenter plane
                // Distance from scoring plane to isocenter along w direction
                double projScale = SAD / (scoringPlaneZ * 10.0); // convert scoringPlaneZ from cm to mm ratio
                double projX = p.position[0] * projScale;
                double projY = p.position[1] * projScale;

                // Check if within jaw opening
                // Jaw defines a symmetric opening: -jawX1 to +jawX2 in X, -jawY1 to +jawY2 in Y
                return (projX < -jawX1 || projX > jawX2 ||
                        projY < -jawY1 || projY > jawY2);
            }),
        m_particles.end());
}

void PhaseSpaceData::filterByType(ParticleType type) {
    m_particles.erase(
        std::remove_if(m_particles.begin(), m_particles.end(),
            [type](const Particle& p) { return p.type != type; }),
        m_particles.end());
}

PhaseSpaceMetrics PhaseSpaceData::computeMetrics() const {
    PhaseSpaceMetrics m;
    if (m_particles.empty()) return m;

    m.totalCount = static_cast<int64_t>(m_particles.size());
    m.minEnergy = m_particles[0].energy;
    m.maxEnergy = m_particles[0].energy;
    m.xRange = {m_particles[0].position[0], m_particles[0].position[0]};
    m.yRange = {m_particles[0].position[1], m_particles[0].position[1]};
    m.zRange = {m_particles[0].position[2], m_particles[0].position[2]};

    double sumE = 0.0, sumU = 0.0, sumV = 0.0;
    double sumU2 = 0.0, sumV2 = 0.0;

    for (const auto& p : m_particles) {
        switch (p.type) {
            case ParticleType::Photon:   ++m.photonCount;   break;
            case ParticleType::Electron: ++m.electronCount; break;
            case ParticleType::Positron: ++m.positronCount; break;
        }

        sumE += p.energy;
        m.minEnergy = std::min(m.minEnergy, p.energy);
        m.maxEnergy = std::max(m.maxEnergy, p.energy);

        sumU += p.direction[0];
        sumV += p.direction[1];
        sumU2 += p.direction[0] * p.direction[0];
        sumV2 += p.direction[1] * p.direction[1];

        m.xRange[0] = std::min(m.xRange[0], p.position[0]);
        m.xRange[1] = std::max(m.xRange[1], p.position[0]);
        m.yRange[0] = std::min(m.yRange[0], p.position[1]);
        m.yRange[1] = std::max(m.yRange[1], p.position[1]);
        m.zRange[0] = std::min(m.zRange[0], p.position[2]);
        m.zRange[1] = std::max(m.zRange[1], p.position[2]);
    }

    double n = static_cast<double>(m.totalCount);
    m.meanEnergy = sumE / n;
    double meanU = sumU / n;
    double meanV = sumV / n;
    m.angularSpreadU = std::sqrt(std::max(0.0, sumU2 / n - meanU * meanU));
    m.angularSpreadV = std::sqrt(std::max(0.0, sumV2 / n - meanV * meanV));

    return m;
}

int64_t PhaseSpaceData::countByType(ParticleType type) const {
    return std::count_if(m_particles.begin(), m_particles.end(),
                         [type](const Particle& p) { return p.type == type; });
}

} // namespace optirad
