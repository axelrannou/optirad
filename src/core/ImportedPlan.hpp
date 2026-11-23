#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace optirad {

/// One DICOM RTPLAN control point, with values carried forward from earlier control points when absent.
/// Positions are in IEC BEAM LIMITING DEVICE coordinates at the isocenter plane (mm).
struct ImportedControlPoint {
    double gantryAngle = 0.0;        // degrees, [0, 360)
    double cumulativeWeight = 0.0;   // normalised to [0, 1] by the beam's FinalCumulativeMetersetWeight
    std::vector<double> bankA;       // MLCX bank A (X1 side) tip positions, one per leaf pair
    std::vector<double> bankB;       // MLCX bank B (X2 side) tip positions
    std::array<double, 2> jawX = {-1e4, 1e4};
    std::array<double, 2> jawY = {-1e4, 1e4};
};

struct ImportedBeam {
    int number = 0;
    std::string name;
    std::string machineName;
    std::string radiationType;       // PHOTON, ELECTRON, ...
    bool isArc = false;
    bool clockwise = true;           // GantryRotationDirection CW (gantry angle increases)
    double sad = 1000.0;             // mm
    double nominalEnergy = 0.0;      // MV
    bool fff = false;                // PrimaryFluenceMode NON_STANDARD / FFF
    double doseRateMUPerMin = 0.0;   // DoseRateSet of the first control point (0 if absent)
    double collimatorAngle = 0.0;    // degrees, first control-point value (dynamic rotation is warned about)
    double couchAngle = 0.0;         // degrees, IEC PatientSupportAngle normalised to (-180, 180]
    std::array<double, 3> isocenter = {0.0, 0.0, 0.0}; // patient LPS mm
    std::vector<double> leafBoundaries; // MLCX LeafPositionBoundaries (pairs + 1)
    double metersetMU = 0.0;         // BeamMeterset: machine MU per fraction
    std::vector<ImportedControlPoint> controlPoints;
};

/// Treatment parameters read from a DICOM RTPLAN (setup beams excluded).
struct ImportedPlan {
    std::string label;
    std::string name;
    std::string sopInstanceUID;
    std::string patientPosition = "HFS";
    int numFractions = 0;
    double prescribedDose = 0.0;     // Gy, whole course (0 if the plan has no target prescription)
    std::vector<ImportedBeam> beams;
    std::vector<std::string> warnings;

    bool hasArcs() const {
        for (const auto& b : beams) if (b.isArc) return true;
        return false;
    }

    /// Isocenter of the first treatment beam (plans with several isocenters are reported by the converter).
    std::optional<std::array<double, 3>> isocenter() const {
        if (beams.empty()) return std::nullopt;
        return beams.front().isocenter;
    }

    bool hasMultipleIsocenters(double tolMm = 0.5) const {
        for (const auto& b : beams)
            for (int i = 0; i < 3; ++i)
                if (std::abs(b.isocenter[i] - beams.front().isocenter[i]) > tolMm) return true;
        return false;
    }
};

} // namespace optirad
