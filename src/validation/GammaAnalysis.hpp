#pragma once

#include "dose/DoseMatrix.hpp"
#include <vector>

namespace optirad {

struct GammaOptions {
    double doseDifferencePercent = 3.0;     // dose criterion [%]
    double distanceToAgreementMm = 2.0;     // DTA criterion [mm]
    bool   globalNormalization   = true;    // global: % of reference max; local: % of local reference dose
    double lowDoseCutoffPercent  = 10.0;    // reference voxels below this % of reference max are excluded
    double searchRadiusFactor    = 2.0;     // search radius = factor * DTA; gamma above factor is reported as factor
    double searchStepFraction    = 0.1;     // search step = fraction * DTA (sub-voxel, trilinear interpolation)
};

struct GammaResult {
    /// Gamma per reference voxel (flat, reference grid order). NaN for voxels excluded by the cutoff.
    std::vector<double> gamma;

    size_t numEvaluated = 0;
    size_t numPassed    = 0;
    double passRate     = 0.0;   // [%]
    double meanGamma    = 0.0;
    double maxGamma     = 0.0;
    double gamma95      = 0.0;   // 95th percentile of evaluated gamma values
};

/// 3D gamma index (Low et al. 1998) of `evaluated` against `reference`.
/// The grids may differ; evaluated dose is trilinearly interpolated at patient-space positions.
class GammaAnalysis {
public:
    static GammaResult compute(const DoseMatrix& reference,
                               const DoseMatrix& evaluated,
                               const GammaOptions& options = {});
};

} // namespace optirad
