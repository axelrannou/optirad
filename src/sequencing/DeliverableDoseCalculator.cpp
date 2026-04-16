#include "DeliverableDoseCalculator.hpp"
#include "core/FluenceMap.hpp"
#include <cmath>
#include <algorithm>

namespace optirad {

DeliverableDoseResult DeliverableDoseCalculator::compute(
    const std::vector<LeafSequenceResult>& beamSequences,
    const Stf& stf,
    const DoseInfluenceMatrix& dij,
    const MachineGeometry& mlc,
    const Grid& doseGrid) {

    DeliverableDoseResult result;
    size_t totalBixels = dij.getNumBixels();
    result.deliverableWeights.assign(totalBixels, 0.0);

    // Compute global bixel offsets per beam (same indexing as PencilBeamEngine)
    size_t numBeams = stf.getCount();
    std::vector<size_t> beamBixelOffset(numBeams + 1, 0);
    for (size_t bi = 0; bi < numBeams; ++bi) {
        const auto* beam = stf.getBeam(bi);
        beamBixelOffset[bi + 1] = beamBixelOffset[bi] + beam->getTotalNumOfBixels();
    }

    // Reconstruct deliverable weights from quantized fluence maps
    // (matRad approach: w = D_0(bixelIndex) / numLevels * calFac)
    for (size_t bi = 0; bi < numBeams && bi < beamSequences.size(); ++bi) {
        const auto& seqResult = beamSequences[bi];
        const auto* beam = stf.getBeam(bi);
        if (!beam || seqResult.quantizedFluence.empty()) continue;

        size_t offset = beamBixelOffset[bi];
        double bw = beam->getBixelWidth();
        double calFac = seqResult.calFac;
        int numLevels = seqResult.numLevels;
        if (numLevels <= 0 || calFac <= 0.0) continue;

        // Reconstruct the fluence map grid parameters from ray positions
        double xMin = 1e30, zMin = 1e30;
        double xMax = -1e30, zMax = -1e30;
        for (size_t ri = 0; ri < beam->getNumOfRays(); ++ri) {
            const auto& pos = beam->getRay(ri)->getRayPosBev();
            xMin = std::min(xMin, pos[0]);
            xMax = std::max(xMax, pos[0]);
            zMin = std::min(zMin, pos[2]);
            zMax = std::max(zMax, pos[2]);
        }
        int nCols = static_cast<int>(std::round((xMax - xMin) / bw)) + 1;
        int nRows = static_cast<int>(std::round((zMax - zMin) / bw)) + 1;

        // For each ray, look up its quantized fluence and convert to weight
        for (size_t ri = 0; ri < beam->getNumOfRays(); ++ri) {
            const auto& pos = beam->getRay(ri)->getRayPosBev();
            int col = static_cast<int>(std::round((pos[0] - xMin) / bw));
            int row = static_cast<int>(std::round((pos[2] - zMin) / bw));

            if (row < 0 || row >= nRows || col < 0 || col >= nCols) continue;

            size_t fmIdx = static_cast<size_t>(row) * nCols + col;
            if (fmIdx >= seqResult.quantizedFluence.size()) continue;

            int D0 = seqResult.quantizedFluence[fmIdx];
            double deliverableWeight = static_cast<double>(D0) / numLevels * calFac;

            size_t globalIdx = offset + ri;
            if (globalIdx < totalBixels) {
                result.deliverableWeights[globalIdx] = deliverableWeight;
            }
        }
    }

    // Compute dose = Dij x deliverableWeights
    auto doseData = dij.computeDose(result.deliverableWeights);
    auto doseMat = std::make_shared<DoseMatrix>();
    doseMat->setGrid(doseGrid);
    doseMat->allocate();
    std::copy(doseData.begin(), doseData.end(), doseMat->data());
    result.dose = std::move(doseMat);

    return result;
}

} // namespace optirad
