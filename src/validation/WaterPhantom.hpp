#pragma once

#include "core/PatientData.hpp"
#include "dose/DoseMatrix.hpp"
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace optirad {

/// Homogeneous water phantom with a single open square field from gantry 0 (beam along +y).
/// Coordinates: isocenter at the origin, water surface at y = -isocenterDepthMm.
struct WaterPhantomConfig {
    double resolutionMm       = 2.5;    // CT and dose grid spacing
    double fieldSizeMm        = 100.0;  // square field at isocenter
    double bixelWidthMm       = 5.0;
    double isocenterDepthMm   = 100.0;  // SSD = SAD - isocenterDepthMm
    double lateralHalfSizeMm  = 200.0;
    double depthMm            = 400.0;
    std::string machineName   = "Generic";
    double relativeThreshold  = 0.0;    // Dij relative threshold (0 = keep everything, as matRad)

    /// Bixel weight as a function of the bixel centre (x, z) at the isocenter plane [mm].
    /// One dose is computed per entry from the same Dij. Empty = uniform weights of 1.
    std::vector<std::function<double(double, double)>> fluences;
};

struct WaterPhantomResult {
    std::shared_ptr<DoseMatrix> dose;                  // first entry of `doses`
    std::vector<std::shared_ptr<DoseMatrix>> doses;    // one per configured fluence
    double surfaceY = 0.0;                             // y of the water surface [mm]
    size_t dijNonZeros = 0;
    double dijSeconds = 0.0;
};

class WaterPhantom {
public:
    /// Throws on failure.
    static WaterPhantomResult compute(const WaterPhantomConfig& config = {});
};

/// A sampled 1D curve: position along the sampled axis [mm] and dose.
struct DoseCurve {
    std::vector<double> position;
    std::vector<double> dose;
};

/// Dose sampled along y (depth) at fixed lateral (x, z), from yStart to yEnd with the given step.
DoseCurve sampleDepthDose(const DoseMatrix& dose, double x, double z,
                          double yStart, double yEnd, double stepMm);

/// Dose sampled along x at fixed (y, z), from xStart to xEnd with the given step.
DoseCurve sampleProfileX(const DoseMatrix& dose, double y, double z,
                         double xStart, double xEnd, double stepMm);

} // namespace optirad
