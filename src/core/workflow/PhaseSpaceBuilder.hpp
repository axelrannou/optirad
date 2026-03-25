#pragma once

#include "core/Plan.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include <memory>
#include <vector>
#include <array>

namespace optirad {

/// Options for loading phase-space beam sources.
struct PhaseSpaceBuildOptions {
    double collimatorAngle = 0.0;
    int64_t maxParticles = 1000000;
    int64_t vizSampleSize = 100000;
};

/// Extracts the shared PSF beam loading logic from CLI loadPhaseSpace()
/// and GUI PhaseSpacePanel: configure → build all beams (OpenMP parallel).
class PhaseSpaceBuilder {
public:
    /// Build phase-space beam sources for all gantry angles in the plan.
    /// Returns one PhaseSpaceBeamSource per beam. Throws on error.
    static std::vector<std::shared_ptr<PhaseSpaceBeamSource>> build(
        const Plan& plan,
        const PhaseSpaceBuildOptions& options = {});
};

} // namespace optirad
