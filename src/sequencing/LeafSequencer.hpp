#pragma once

#include "core/Aperture.hpp"
#include "core/FluenceMap.hpp"
#include "core/Machine.hpp"
#include <vector>

namespace optirad {

/// Step-and-shoot leaf sequencer: converts a 2D fluence map into a sequence
/// of MLC apertures using intensity level decomposition.
class LeafSequencer {
public:
    /// Decompose a fluence map into step-and-shoot aperture segments.
    /// @param fluence       2D fluence map for one beam
    /// @param mlc           Machine MLC geometry (leaf widths, max travel, interdigitation)
    /// @param opts          Sequencing options (num levels, min segment MU, etc.)
    /// @return              Sequencing result with apertures and quality metrics
    static LeafSequenceResult sequenceBeam(
        const FluenceMap& fluence,
        const MachineGeometry& mlc,
        const LeafSequencerOptions& opts);
};

} // namespace optirad
