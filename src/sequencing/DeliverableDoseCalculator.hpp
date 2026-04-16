#pragma once

#include "core/Aperture.hpp"
#include "core/FluenceMap.hpp"
#include "core/Stf.hpp"
#include "core/Machine.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include <vector>
#include <memory>

namespace optirad {

/// Result of deliverable dose computation.
struct DeliverableDoseResult {
    /// Effective bixel weights after leaf sequencing.
    std::vector<double> deliverableWeights;

    /// Dose matrix computed from deliverable weights.
    std::shared_ptr<DoseMatrix> dose;
};

/// Converts leaf sequence apertures back into effective bixel weights
/// and computes the resulting deliverable dose via the Dij matrix.
class DeliverableDoseCalculator {
public:
    /// Compute deliverable dose from leaf sequences.
    /// For each beam, determines which bixels are open in each segment,
    /// builds effective weights, and computes dose = Dij x weights.
    static DeliverableDoseResult compute(
        const std::vector<LeafSequenceResult>& beamSequences,
        const Stf& stf,
        const DoseInfluenceMatrix& dij,
        const MachineGeometry& mlc,
        const Grid& doseGrid);
};

} // namespace optirad
