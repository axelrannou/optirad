#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace optirad {

/// Particle type classification from IAEA phase-space data
enum class ParticleType : int8_t {
    Photon = 1,
    Electron = -1,
    Positron = 2
};

/// A single particle record from an IAEA phase-space file.
struct Particle {
    ParticleType type = ParticleType::Photon;
    std::array<double, 3> position = {0.0, 0.0, 0.0};   // (x, y, z) in mm
    std::array<double, 3> direction = {0.0, 0.0, 1.0};   // (u, v, w) direction cosines
    double energy = 0.0;   // MeV
    double weight = 1.0;
};

/// Verification metrics for a phase-space dataset
struct PhaseSpaceMetrics {
    int64_t totalCount = 0;
    int64_t photonCount = 0;
    int64_t electronCount = 0;
    int64_t positronCount = 0;

    double minEnergy = 0.0;
    double maxEnergy = 0.0;
    double meanEnergy = 0.0;

    double angularSpreadU = 0.0; // std dev of u direction cosine
    double angularSpreadV = 0.0; // std dev of v direction cosine

    // Spatial extent (mm)
    std::array<double, 2> xRange = {0.0, 0.0};
    std::array<double, 2> yRange = {0.0, 0.0};
    std::array<double, 2> zRange = {0.0, 0.0};
};

/// Container for phase-space particle data with filtering, sampling, and statistics.
class PhaseSpaceData {
public:
    PhaseSpaceData() = default;

    /// Direct access to the particle vector
    std::vector<Particle>& particles() { return m_particles; }
    const std::vector<Particle>& particles() const { return m_particles; }

    /// Add a single particle
    void addParticle(const Particle& p) { m_particles.push_back(p); }

    /// Reserve space for expected number of particles
    void reserve(size_t n) { m_particles.reserve(n); }

    /// Number of particles
    size_t size() const { return m_particles.size(); }
    bool empty() const { return m_particles.empty(); }

    /// Clear all particles
    void clear() { m_particles.clear(); }

    /// Return a random subset of n particles for visualization
    PhaseSpaceData sample(size_t n, unsigned int seed = 42) const;

    /// Remove particles outside the jaw aperture.
    /// Jaw positions define a rectangular opening at the isocenter plane.
    /// Particles are projected to the isocenter plane using direction cosines.
    /// @param jawX1, jawX2: X jaw positions (mm, positive = open from center)
    /// @param jawY1, jawY2: Y jaw positions (mm, positive = open from center)
    /// @param SAD: source-axis distance (mm) for projection
    /// @param scoringPlaneZ: Z position of scoring plane (mm)
    void filterByJaws(double jawX1, double jawX2, double jawY1, double jawY2,
                      double SAD, double scoringPlaneZ);

    /// Filter by particle type (keep only the specified type)
    void filterByType(ParticleType type);

    /// Compute verification metrics for the current particle set
    PhaseSpaceMetrics computeMetrics() const;

    /// Get count by particle type
    int64_t countByType(ParticleType type) const;

private:
    std::vector<Particle> m_particles;
};

} // namespace optirad
