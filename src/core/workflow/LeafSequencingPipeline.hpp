#pragma once

#include "core/Aperture.hpp"
#include "core/FluenceMap.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/PatientData.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "core/workflow/PlanAnalysis.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <vector>
#include <string>
#include <functional>

namespace optirad {

/// Result of the leaf sequencing pipeline.
struct LeafSequencingPipelineResult {
    /// Per-beam sequencing results (apertures, MU, fidelity).
    std::vector<LeafSequenceResult> beamSequences;

    /// Effective bixel weights after aperture decomposition.
    std::vector<double> deliverableWeights;

    /// Deliverable dose matrix.
    std::shared_ptr<DoseMatrix> deliverableDose;

    /// Per-structure statistics for the deliverable dose.
    std::vector<StructureDoseStats> deliverableStats;

    /// Total segments across all beams.
    int totalSegments = 0;

    /// Total MU across all beams.
    double totalMU = 0.0;

    /// Mean fluence fidelity across beams (Pearson correlation).
    double meanFidelity = 0.0;
};

/// Orchestrates the leaf sequencing workflow:
///   optimized weights → fluence maps → sequencing → deliverable dose → stats.
class LeafSequencingPipeline {
public:
    using ProgressCallback = std::function<void(int current, int total)>;

    /// Run the leaf sequencing pipeline.
    /// @param weights       Optimized bixel weights from the optimizer
    /// @param stf           Steering file with beam geometry
    /// @param dij           Dose influence matrix
    /// @param plan          Treatment plan (for machine geometry)
    /// @param patientData   Patient data (for structure-based stats)
    /// @param doseGrid      Dose calculation grid
    /// @param opts          Leaf sequencer options
    /// @param progressCb    Optional progress callback (beamIndex, totalBeams)
    static LeafSequencingPipelineResult run(
        const std::vector<double>& weights,
        const Stf& stf,
        const DoseInfluenceMatrix& dij,
        const Plan& plan,
        const PatientData& patientData,
        const Grid& doseGrid,
        const LeafSequencerOptions& opts,
        ProgressCallback progressCb = nullptr);
};

} // namespace optirad
