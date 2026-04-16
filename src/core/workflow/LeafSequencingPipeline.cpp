#include "LeafSequencingPipeline.hpp"
#include "sequencing/LeafSequencer.hpp"
#include "sequencing/DeliverableDoseCalculator.hpp"
#include "utils/Logger.hpp"
#include <chrono>
#include <numeric>

namespace optirad {

LeafSequencingPipelineResult LeafSequencingPipeline::run(
    const std::vector<double>& weights,
    const Stf& stf,
    const DoseInfluenceMatrix& dij,
    const Plan& plan,
    const PatientData& patientData,
    const Grid& doseGrid,
    const LeafSequencerOptions& opts,
    ProgressCallback progressCb) {

    LeafSequencingPipelineResult result;
    auto startTime = std::chrono::steady_clock::now();

    const auto& mlc = plan.getMachine().getGeometry();
    size_t numBeams = stf.getCount();

    Logger::info("LeafSequencingPipeline: " + std::to_string(numBeams) +
                 " beams, " + std::to_string(opts.numLevels) + " intensity levels");

    // Compute global bixel offsets
    std::vector<size_t> beamBixelOffset(numBeams + 1, 0);
    for (size_t bi = 0; bi < numBeams; ++bi) {
        beamBixelOffset[bi + 1] = beamBixelOffset[bi] +
                                   stf.getBeam(bi)->getTotalNumOfBixels();
    }

    // Sequence each beam
    result.beamSequences.resize(numBeams);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t bi = 0; bi < numBeams; ++bi) {
        const auto* beam = stf.getBeam(bi);
        auto fluence = FluenceMap::fromBeamWeights(*beam, weights, beamBixelOffset[bi]);
        result.beamSequences[bi] = LeafSequencer::sequenceBeam(fluence, mlc, opts);
        result.beamSequences[bi].beamIndex = bi;

        if (progressCb) {
            progressCb(static_cast<int>(bi + 1), static_cast<int>(numBeams));
        }
    }

    // Aggregate statistics
    result.totalSegments = 0;
    result.totalMU = 0.0;
    double fidelitySum = 0.0;
    int fidelityCount = 0;
    for (const auto& seq : result.beamSequences) {
        result.totalSegments += static_cast<int>(seq.segments.size());
        result.totalMU += seq.totalMU;
        if (!seq.segments.empty()) {
            fidelitySum += seq.fluenceFidelity;
            ++fidelityCount;
        }
    }
    result.meanFidelity = (fidelityCount > 0) ? fidelitySum / fidelityCount : 0.0;

    Logger::info("LeafSequencingPipeline: " + std::to_string(result.totalSegments) +
                 " segments, " + std::to_string(result.totalMU) + " MU, fidelity=" +
                 std::to_string(result.meanFidelity));

    // Compute deliverable dose
    auto delResult = DeliverableDoseCalculator::compute(
        result.beamSequences, stf, dij, mlc, doseGrid);
    result.deliverableWeights = std::move(delResult.deliverableWeights);
    result.deliverableDose = std::move(delResult.dose);

    // Compute plan statistics for the deliverable dose
    double prescribedDose = plan.getPrescribedDose();
    if (prescribedDose <= 0.0) prescribedDose = 66.0; // fallback if not set

    if (result.deliverableDose) {
        result.deliverableStats = PlanAnalysis::computeStats(
            *result.deliverableDose, patientData, doseGrid, prescribedDose);

        Logger::info("LeafSequencingPipeline: deliverable dose max=" +
                     std::to_string(result.deliverableDose->getMax()) + " Gy");
    }

    auto endTime = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(endTime - startTime).count();
    Logger::info("LeafSequencingPipeline: completed in " + std::to_string(elapsed) + "s");

    return result;
}

} // namespace optirad
