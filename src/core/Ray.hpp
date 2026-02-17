#pragma once

#include "../geometry/MathUtils.hpp"
#include <array>
#include <vector>

namespace optirad {

/**
 * Ray - represents a single beamlet/ray within a beam.
 * 
 * Each ray has a position and target point in both beam's eye view (BEV)
 * and the patient (LPS) coordinate system. For photon beams, rays also 
 * store the 4-corner geometry at the isocenter plane and at the 
 * source-to-collimator distance (SCD).
 */
class Ray {
public:
    Ray() = default;

    // --- BEV (Beam's Eye View) coordinates ---
    void setRayPosBev(const Vec3& pos) { m_rayPosBev = pos; }
    const Vec3& getRayPosBev() const { return m_rayPosBev; }

    void setTargetPointBev(const Vec3& target) { m_targetPointBev = target; }
    const Vec3& getTargetPointBev() const { return m_targetPointBev; }

    // --- LPS (patient) coordinates ---
    void setRayPos(const Vec3& pos) { m_rayPos = pos; }
    const Vec3& getRayPos() const { return m_rayPos; }

    void setTargetPoint(const Vec3& target) { m_targetPoint = target; }
    const Vec3& getTargetPoint() const { return m_targetPoint; }

    // --- Photon-specific: beamlet corners at isocenter plane (4 corners, 3D each) ---
    void setBeamletCornersAtIso(const std::array<Vec3, 4>& corners) { m_beamletCornersAtIso = corners; }
    const std::array<Vec3, 4>& getBeamletCornersAtIso() const { return m_beamletCornersAtIso; }

    // --- Photon-specific: ray corners projected to SCD plane (4 corners, 3D each) ---
    void setRayCornersSCD(const std::array<Vec3, 4>& corners) { m_rayCornersSCD = corners; }
    const std::array<Vec3, 4>& getRayCornersSCD() const { return m_rayCornersSCD; }

    // --- Energy (MV for photons, MeV for particles) ---
    void setEnergy(double energy) { m_energy = energy; }
    double getEnergy() const { return m_energy; }

    // --- Number of bixels for this ray (1 for photons, multiple for particles) ---
    size_t getNumOfBixels() const { return 1; } // Photon IMRT: always 1 bixel per ray

private:
    // BEV coordinates
    Vec3 m_rayPosBev = {0.0, 0.0, 0.0};
    Vec3 m_targetPointBev = {0.0, 0.0, 0.0};

    // LPS (patient) coordinates
    Vec3 m_rayPos = {0.0, 0.0, 0.0};
    Vec3 m_targetPoint = {0.0, 0.0, 0.0};

    // Photon-specific corner geometry
    std::array<Vec3, 4> m_beamletCornersAtIso = {Vec3{0,0,0}, Vec3{0,0,0}, Vec3{0,0,0}, Vec3{0,0,0}};
    std::array<Vec3, 4> m_rayCornersSCD = {Vec3{0,0,0}, Vec3{0,0,0}, Vec3{0,0,0}, Vec3{0,0,0}};

    // Energy
    double m_energy = 6.0; // Default photon energy (MV)
};

} // namespace optirad
