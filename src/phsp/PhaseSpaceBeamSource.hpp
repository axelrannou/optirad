#pragma once

#include "phsp/PhaseSpaceData.hpp"
#include "phsp/IAEAHeaderParser.hpp"
#include "core/Machine.hpp"
#include "geometry/MathUtils.hpp"

#include <string>
#include <vector>
#include <array>
#include <memory>

namespace optirad {

/// Beam source built from IAEA phase-space files.
/// Replaces STF for phase-space machines — holds particle-based beam data
/// rather than ray-based Beam/Stf objects.
///
/// Pipeline: load PSF → filter by aperture → transform to patient coords → sample for viz
class PhaseSpaceBeamSource {
public:
    PhaseSpaceBeamSource() = default;

    /// Configure the beam source from a phase-space machine.
    /// @param machine       Machine with MachineGeometry (PSF paths, jaws, MLC, etc.)
    /// @param gantryAngle   Gantry angle in degrees
    /// @param collimatorAngle Collimator angle in degrees
    /// @param couchAngle    Couch angle in degrees
    /// @param isocenter     Isocenter position in LPS (mm)
    void configure(const Machine& machine,
                   double gantryAngle,
                   double collimatorAngle,
                   double couchAngle,
                   const std::array<double, 3>& isocenter);

    /// Build the beam source: read PSF files, filter, transform.
    /// @param maxParticles  Max particles to read (0 = all from first file)
    /// @param vizSampleSize Number of particles to keep for visualization
    void build(int64_t maxParticles = 0, int64_t vizSampleSize = 100000);

    /// Whether the source has been built successfully
    bool isBuilt() const { return m_built; }

    /// Get the full (filtered, transformed) phase-space data
    const PhaseSpaceData& getData() const { return m_data; }

    /// Get the visualization sample (subset of particles for rendering)
    const PhaseSpaceData& getVisualizationSample() const { return m_vizSample; }

    /// Get verification metrics
    const PhaseSpaceMetrics& getMetrics() const { return m_metrics; }

    /// Get the isocenter position (LPS, mm)
    const std::array<double, 3>& getIsocenter() const { return m_isocenter; }

    /// Get the source position (LPS, mm) — computed from SAD + gantry rotation
    const std::array<double, 3>& getSourcePosition() const { return m_sourcePos; }

    /// Get gantry angle
    double getGantryAngle() const { return m_gantryAngle; }

    /// Get the IAEA header info (available after build)
    const IAEAHeaderInfo& getHeaderInfo() const { return m_headerInfo; }

    /// Compute energy histogram bins.
    /// Returns pairs of (bin_center_MeV, count).
    std::vector<std::pair<double, int64_t>> computeEnergyHistogram(int numBins = 50) const;

private:
    /// Transform particle positions/directions from PSF frame → patient LPS
    void transformToPatientCoords();

    /// Compute source position in LPS from SAD + gantry/couch rotation
    void computeSourcePosition();

    // Configuration
    Machine m_machine;
    double m_gantryAngle = 0.0;
    double m_collimatorAngle = 0.0;
    double m_couchAngle = 0.0;
    std::array<double, 3> m_isocenter = {0.0, 0.0, 0.0};

    // Computed
    std::array<double, 3> m_sourcePos = {0.0, 0.0, 0.0};
    Mat3 m_rotMatrix;

    // Data
    PhaseSpaceData m_data;         // full filtered + transformed data
    PhaseSpaceData m_vizSample;    // subset for visualization
    PhaseSpaceMetrics m_metrics;
    IAEAHeaderInfo m_headerInfo;

    bool m_built = false;
};

} // namespace optirad
