#pragma once

#include "Ray.hpp"
#include "../geometry/MathUtils.hpp"
#include <string>
#include <array>
#include <vector>

namespace optirad {

/**
 * Beam - represents a single treatment beam with its geometry and ray collection.
 *
 * A beam has angular parameters (gantry/couch), geometric parameters (isocenter, SAD, SCD),
 * and contains a collection of rays (beamlets). Each ray defines a pencil beam direction
 * and its geometry in both BEV and patient (LPS) coordinate systems.
 */
class Beam {
public:
    Beam() = default;

    // --- Angular parameters ---
    void setGantryAngle(double angle);
    void setCouchAngle(double angle);
    double getGantryAngle() const;
    double getCouchAngle() const;

    // --- Isocenter ---
    void setIsocenter(double x, double y, double z);
    void setIsocenter(const Vec3& iso);
    const Vec3& getIsocenter() const { return m_isocenter; }

    // --- Bixel width ---
    void setBixelWidth(double width) { m_bixelWidth = width; }
    double getBixelWidth() const { return m_bixelWidth; }

    // --- Radiation mode ---
    void setRadiationMode(const std::string& mode) { m_radiationMode = mode; }
    const std::string& getRadiationMode() const { return m_radiationMode; }

    // --- Machine name ---
    void setMachineName(const std::string& name) { m_machineName = name; }
    const std::string& getMachineName() const { return m_machineName; }

    // --- SAD (Source-Axis Distance) ---
    void setSAD(double sad) { m_SAD = sad; }
    double getSAD() const { return m_SAD; }

    // --- SCD (Source-Collimator Distance) ---
    void setSCD(double scd) { m_SCD = scd; }
    double getSCD() const { return m_SCD; }

    // --- Source point in BEV ---
    const Vec3& getSourcePointBev() const { return m_sourcePointBev; }

    // --- Source point in LPS ---
    const Vec3& getSourcePoint() const { return m_sourcePoint; }

    // --- Ray management ---
    void addRay(const Ray& ray);
    void addRay(Ray&& ray);
    const Ray* getRay(size_t index) const;
    Ray* getRay(size_t index);
    const std::vector<Ray>& getRays() const { return m_rays; }
    std::vector<Ray>& rays() { return m_rays; }
    size_t getNumOfRays() const { return m_rays.size(); }

    // --- Bixel counts ---
    std::vector<size_t> getNumOfBixelsPerRay() const;
    size_t getTotalNumOfBixels() const;

    // --- Source point computation ---
    /// Computes source points in BEV and LPS based on SAD and angles
    void computeSourcePoints();

    // --- Ray initialization ---
    /// Initializes ray positions: given a set of BEV positions, creates rays
    /// with proper target points and transforms to LPS coordinates
    void initRaysFromPositions(const std::vector<Vec3>& rayPositionsBev);

    // --- Photon-specific ray corner computation ---
    /// Computes beamlet corners at isocenter and SCD planes for all rays
    void computePhotonRayCorners();

    // --- Energy assignment ---
    /// Assigns the given energy to all rays
    void setAllRayEnergies(double energy);

    // --- Field size ---
    void setFieldSize(const std::array<double, 2>& fieldSize) { m_fieldSize = fieldSize; }
    const std::array<double, 2>& getFieldSize() const { return m_fieldSize; }

    // --- Simple ray generation (grid-based without CT) ---
    /// Generates rays on a regular bixel grid covering the given field size
    void generateRays(double bixelWidth, const std::array<double, 2>& fieldSize);

    // --- Target-aware ray generation ---
    /// Generates rays by projecting target voxel world coordinates onto the BEV
    /// isocenter plane, snapping to the bixel grid, and keeping unique positions.
    /// This is the proper method that produces angle-dependent ray counts.
    /// The target world coordinates should already include any 3D margin expansion.
    /// @param targetWorldCoords  Target voxel positions in LPS world coordinates
    /// @param bixelWidth         Bixel width in mm
    /// @param ctResolution       CT voxel resolution {dy, dx, dz} in mm for additional padding
    void generateRaysFromTarget(const std::vector<Vec3>& targetWorldCoords,
                                double bixelWidth,
                                const Vec3& ctResolution = {1.0, 1.0, 1.0});

private:
    // Angular parameters
    double m_gantryAngle = 0.0;
    double m_couchAngle = 0.0;

    // Geometry
    Vec3 m_isocenter = {0.0, 0.0, 0.0};
    double m_bixelWidth = 7.0;  // mm
    double m_SAD = 1000.0;      // mm
    double m_SCD = 500.0;       // mm

    // Source points
    Vec3 m_sourcePointBev = {0.0, 0.0, 0.0};
    Vec3 m_sourcePoint = {0.0, 0.0, 0.0};

    // Metadata
    std::string m_radiationMode = "photons";
    std::string m_machineName = "Generic";

    // Field size [width, height] in mm
    std::array<double, 2> m_fieldSize = {100.0, 100.0};

    // Rays
    std::vector<Ray> m_rays;
};

} // namespace optirad
