#include "PlanAnalysis.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

double PlanAnalysis::computeDx(const std::vector<double>& sortedDoses, double percentile) {
    if (sortedDoses.empty()) return 0.0;
    size_t N = sortedDoses.size();
    if (N == 1) return sortedDoses[0];
    // Dx% = dose such that x% of the volume receives at least that dose.
    // Use linear interpolation between ranks (matches MATLAB prctile / numpy percentile).
    double pos = (1.0 - percentile / 100.0) * static_cast<double>(N - 1);
    size_t lo = static_cast<size_t>(pos);
    if (lo >= N - 1) return sortedDoses[N - 1];
    double frac = pos - static_cast<double>(lo);
    return sortedDoses[lo] + frac * (sortedDoses[lo + 1] - sortedDoses[lo]);
}

double PlanAnalysis::computeVx(const std::vector<double>& sortedDoses, double doseThreshold) {
    if (sortedDoses.empty()) return 0.0;
    // VxGy = percentage of volume receiving >= x Gy
    // Use lower_bound on sorted ascending array
    auto it = std::lower_bound(sortedDoses.begin(), sortedDoses.end(), doseThreshold);
    size_t countAbove = sortedDoses.end() - it;
    return 100.0 * static_cast<double>(countAbove) / sortedDoses.size();
}

std::vector<StructureDoseStats> PlanAnalysis::computeStats(
    const DoseMatrix& dose,
    const PatientData& patient,
    const Grid& doseGrid,
    double prescribedDose)
{
    std::vector<StructureDoseStats> results;

    const auto* ss = patient.getStructureSet();
    if (!ss) return results;

    const auto& ctGrid = patient.getGrid();
    const double* doseData = dose.data();
    size_t doseSize = dose.size();

    // Dose grid dimensions for bounds checking
    auto doseDims = doseGrid.getDimensions();
    double dNx = static_cast<double>(doseDims[0]);
    double dNy = static_cast<double>(doseDims[1]);
    double dNz = static_cast<double>(doseDims[2]);

    // CT grid dimensions for index decomposition
    auto ctDims = ctGrid.getDimensions();
    size_t ctNx = ctDims[0], ctNy = ctDims[1];

    // For CI: count total voxels receiving >= prescribed dose
    size_t totalVoxelsAboveRx = 0;
    {
        for (size_t i = 0; i < doseSize; ++i) {
            if (doseData[i] >= prescribedDose) ++totalVoxelsAboveRx;
        }
    }

    // Collect non-empty structure indices for OpenMP parallelization
    std::vector<size_t> structIndices;
    for (size_t si = 0; si < ss->getCount(); ++si) {
        const auto* structure = ss->getStructure(si);
        if (structure && !structure->getVoxelIndices().empty())
            structIndices.push_back(si);
    }

    results.resize(structIndices.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < structIndices.size(); ++idx) {
        size_t si = structIndices[idx];
        const auto* structure = ss->getStructure(si);

        const auto& ctVoxels = structure->getVoxelIndices();

        // Collect dose values via trilinear interpolation at each CT-voxel center.
        // This preserves the fine CT-resolution structure boundary instead of
        // rounding to nearest dose voxel (which dilates the boundary).
        std::vector<double> doseValues;
        doseValues.reserve(ctVoxels.size());
        for (size_t idx : ctVoxels) {
            // Decompose flat CT index to (ix, iy, iz)
            size_t iz = idx / (ctNx * ctNy);
            size_t rem = idx % (ctNx * ctNy);
            size_t iy = rem / ctNx;
            size_t ix = rem % ctNx;

            // CT voxel center → patient LPS → fractional dose voxel
            Vec3 lps = ctGrid.voxelToPatient({
                static_cast<double>(ix),
                static_cast<double>(iy),
                static_cast<double>(iz)});
            Vec3 fijk = doseGrid.patientToVoxel(lps);

            // Skip if outside dose grid (with small margin)
            if (fijk[0] < -0.5 || fijk[0] > dNx - 0.5 ||
                fijk[1] < -0.5 || fijk[1] > dNy - 0.5 ||
                fijk[2] < -0.5 || fijk[2] > dNz - 0.5) {
                continue;
            }

            doseValues.push_back(dose.interpolateAt(fijk[0], fijk[1], fijk[2]));
        }
        if (doseValues.empty()) continue;

        // Sort for percentile calculations
        std::sort(doseValues.begin(), doseValues.end());

        StructureDoseStats stats;
        stats.name = structure->getName();
        stats.type = structure->getType();
        stats.numVoxels = doseValues.size();
        stats.sortedDoses = doseValues;

        // Basic stats
        stats.minDose = doseValues.front();
        stats.maxDose = doseValues.back();
        stats.meanDose = std::accumulate(doseValues.begin(), doseValues.end(), 0.0) / doseValues.size();

        // Standard deviation
        double sqSum = 0.0;
        for (double d : doseValues) {
            double diff = d - stats.meanDose;
            sqSum += diff * diff;
        }
        stats.stdDose = std::sqrt(sqSum / doseValues.size());

        // Dx% values
        stats.d2 = computeDx(doseValues, 2.0);
        stats.d5 = computeDx(doseValues, 5.0);
        stats.d50 = computeDx(doseValues, 50.0);
        stats.d95 = computeDx(doseValues, 95.0);
        stats.d98 = computeDx(doseValues, 98.0);

        // VxGy values
        stats.v20 = computeVx(doseValues, 20.0);
        stats.v40 = computeVx(doseValues, 40.0);
        stats.v50 = computeVx(doseValues, 50.0);
        stats.v60 = computeVx(doseValues, 60.0);

        // Quality indices for targets
        bool isTarget = (stats.type.find("TARGET") != std::string::npos ||
                        stats.name.find("PTV") != std::string::npos ||
                        stats.name.find("CTV") != std::string::npos ||
                        stats.name.find("GTV") != std::string::npos);

        if (isTarget) {
            // CI = V_Rx_in_target / V_target (ideally 1.0)
            size_t voxelsAboveRxInTarget = 0;
            for (double d : doseValues) {
                if (d >= prescribedDose) ++voxelsAboveRxInTarget;
            }
            stats.conformityIndex = static_cast<double>(voxelsAboveRxInTarget) / doseValues.size();

            // HI = (D2 - D98) / D50 (ideally 0.0 — homogeneous)
            if (stats.d50 > 0.0) {
                stats.homogeneityIndex = (stats.d2 - stats.d98) / stats.d50;
            }
        }

        results[idx] = std::move(stats);
    }

    return results;
}

std::vector<DVHCurveData> PlanAnalysis::computeDVHCurves(
    const std::vector<StructureDoseStats>& stats,
    double maxDose,
    int numBins)
{
    std::vector<DVHCurveData> curves;

    for (const auto& s : stats) {
        DVHCurveData curve;
        curve.structureName = s.name;
        curve.doses.resize(numBins);
        curve.volumes.resize(numBins);

        double step = maxDose / (numBins - 1);
        const auto& sorted = s.sortedDoses;

        for (int i = 0; i < numBins; ++i) {
            double doseThreshold = i * step;
            curve.doses[i] = static_cast<float>(doseThreshold);

            // Cumulative DVH: % volume receiving >= threshold
            auto it = std::lower_bound(sorted.begin(), sorted.end(), doseThreshold);
            size_t countAbove = sorted.end() - it;
            curve.volumes[i] = static_cast<float>(100.0 * countAbove / sorted.size());
        }

        curves.push_back(std::move(curve));
    }

    return curves;
}

void PlanAnalysis::print(const std::vector<StructureDoseStats>& stats, std::ostream& os) {
    os << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
    os << "║                           PLAN ANALYSIS                                      ║\n";
    os << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";

    for (const auto& s : stats) {
        bool isTarget = (s.type.find("TARGET") != std::string::npos ||
                        s.name.find("PTV") != std::string::npos ||
                        s.name.find("CTV") != std::string::npos ||
                        s.name.find("GTV") != std::string::npos);

        os << "┌─ " << s.name << " (" << s.type << ") ── " << s.numVoxels << " voxels ─┐\n";
        os << "│  Dose:  min = " << std::fixed << std::setprecision(2)
           << std::setw(7) << s.minDose << " Gy  |  max = "
           << std::setw(7) << s.maxDose << " Gy  |  mean = "
           << std::setw(7) << s.meanDose << " Gy  |  std = "
           << std::setw(7) << s.stdDose << " Gy\n";

        os << "│  Dx%:   D2  = " << std::setw(7) << s.d2
           << " Gy  |  D5  = " << std::setw(7) << s.d5
           << " Gy  |  D50 = " << std::setw(7) << s.d50
           << " Gy  |  D95 = " << std::setw(7) << s.d95
           << " Gy  |  D98 = " << std::setw(7) << s.d98 << " Gy\n";

        os << "│  VxGy:  V20 = " << std::setw(6) << std::setprecision(1) << s.v20
           << "%  |  V40 = " << std::setw(6) << s.v40
           << "%  |  V50 = " << std::setw(6) << s.v50
           << "%  |  V60 = " << std::setw(6) << s.v60 << "%\n";

        if (isTarget) {
            os << "│  CI = " << std::setprecision(3) << s.conformityIndex
               << "  |  HI = " << std::setprecision(3) << s.homogeneityIndex << "\n";
        }
        os << "└──────────────────────────────────────────────────────────────────────────────┘\n\n";
    }
}

} // namespace optirad
