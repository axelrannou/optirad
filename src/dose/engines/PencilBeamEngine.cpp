#include "PencilBeamEngine.hpp"
#include "../SiddonRayTracer.hpp"
#include "../SSDCalculator.hpp"
#include "../RadDepthCalculator.hpp"
#include "utils/Logger.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <unordered_set>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

std::string PencilBeamEngine::getName() const {
    return "PencilBeam";
}

// ────────────────────────────────────────────────
// initDoseCalc: cache machine data and set up parameters
// ────────────────────────────────────────────────
void PencilBeamEngine::initDoseCalc(const Plan& plan, const Grid& doseGrid) {
    const auto& machine = plan.getMachine();
    const auto& data = machine.getData();
    const auto& meta = machine.getMeta();

    m_betas = data.betas;
    m_mu = data.m;
    m_SAD = meta.SAD;
    m_SCD = meta.SCD;
    m_penumbraFWHM = data.penumbraFWHMatIso;
    m_kernelPos = data.kernelPos;
    m_kernelEntries = data.kernel;

    // Store full kernel extent (for convolution grid)
    if (!m_kernelPos.empty()) {
        m_kernelCutOff = m_kernelPos.back();
    } else {
        m_kernelCutOff = 50.0;
    }

    // Effective lateral cutoff: geometric cutoff + field width diagonal contribution
    // Matches matRad: effectiveLateralCutOff = geometricLateralCutOff + fieldWidth/sqrt(2)
    m_lateralCutOff = m_geometricLateralCutOff + m_bixelWidth / std::sqrt(2.0);

    Logger::info("PencilBeamEngine init: betas=[" +
        std::to_string(m_betas[0]) + "," + std::to_string(m_betas[1]) + "," + std::to_string(m_betas[2]) +
        "], mu=" + std::to_string(m_mu) +
        ", SAD=" + std::to_string(m_SAD) +
        ", penumbraFWHM=" + std::to_string(m_penumbraFWHM) +
        ", kernelEntries=" + std::to_string(m_kernelEntries.size()) +
        ", kernelCutOff=" + std::to_string(m_kernelCutOff) +
        ", geometricLateralCutOff=" + std::to_string(m_geometricLateralCutOff) +
        ", effectiveLateralCutOff=" + std::to_string(m_lateralCutOff));
}

// ────────────────────────────────────────────────
// initBeam: per-beam setup (rotate voxels to BEV, compute rad depths, SSD, kernel interpolators)
// ────────────────────────────────────────────────
PencilBeamEngine::BeamData PencilBeamEngine::initBeam(
    const Beam& beam,
    const PatientData& patientData,
    const Grid& doseGrid,
    const std::vector<size_t>& allVoxelIndices)
{
    BeamData bd;

    auto dims = doseGrid.getDimensions();
    size_t nx = dims[0], ny = dims[1];

    // Get rotation matrix: LPS → BEV
    Mat3 rotMat = getRotationMatrix(beam.getGantryAngle(), beam.getCouchAngle());
    Mat3 rotMatInv = transpose(rotMat); // BEV → LPS (rotation matrices are orthogonal)

    Vec3 iso = beam.getIsocenter();
    Vec3 sourcePoint = beam.getSourcePoint();

    // Electron density data
    const double* edData = nullptr;
    if (patientData.getEDVolume()) {
        edData = patientData.getEDVolume()->data();
    }
    const auto& ctGrid = patientData.getGrid();

    // 1. Rotate all voxels to BEV and compute geometric distances
    const size_t nVox = allVoxelIndices.size();
    bd.bevCoords.resize(nVox);
    bd.geoDistances.resize(nVox);
    bd.voxelIndices.resize(nVox);

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < nVox; ++i) {
        size_t idx = allVoxelIndices[i];

        // Convert flat index to ijk
        size_t iz = idx / (nx * ny);
        size_t rem = idx % (nx * ny);
        size_t iy = rem / nx;
        size_t ix = rem % nx;

        Vec3 lps = doseGrid.voxelToPatient({
            static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)});

        // Translate to isocenter-centered coordinates then rotate
        Vec3 centered = vecSub(lps, iso);
        Vec3 bev = rotMatInv * centered;

        // Geometric distance from source to voxel
        Vec3 srcToVox = vecSub(lps, sourcePoint);
        double geoDist = norm(srcToVox);

        bd.voxelIndices[i] = idx;
        bd.bevCoords[i] = bev;
        bd.geoDistances[i] = geoDist;
    }

    // 1b. Beam-level pre-filter: only keep voxels within reach of ANY ray in this beam.
    // The beam field spans from the outermost ray positions, so the effective beam-level
    // cutoff is: maxRayOffset + lateralCutOff (projected at isocenter plane).
    // This matches matRad's implicit filtering from matRad_rayTracing coverage.
    {
        // Compute maximum ray offset from beam axis at isocenter plane
        double maxRayOffset = 0.0;
        for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
            const auto* ray = beam.getRay(r);
            if (!ray) continue;
            Vec3 rp = ray->getRayPosBev();
            double rayOff = std::sqrt(rp[0] * rp[0] + rp[2] * rp[2]);
            if (rayOff > maxRayOffset) maxRayOffset = rayOff;
        }
        double beamCutoff = maxRayOffset + m_lateralCutOff;
        double cutoffRatioSq = (beamCutoff / m_SAD) * (beamCutoff / m_SAD);

        size_t writeIdx = 0;
        for (size_t i = 0; i < nVox; ++i) {
            const Vec3& bev = bd.bevCoords[i];
            double depth = m_SAD + bev[1]; // distance from source along beam direction
            if (depth > 1.0) {
                double latDistSq = bev[0] * bev[0] + bev[2] * bev[2];
                if (latDistSq <= cutoffRatioSq * depth * depth) {
                    if (writeIdx != i) {
                        bd.bevCoords[writeIdx] = bd.bevCoords[i];
                        bd.geoDistances[writeIdx] = bd.geoDistances[i];
                        bd.voxelIndices[writeIdx] = bd.voxelIndices[i];
                    }
                    writeIdx++;
                }
            }
        }
        bd.bevCoords.resize(writeIdx);
        bd.geoDistances.resize(writeIdx);
        bd.voxelIndices.resize(writeIdx);
    }

    // 2. Compute radiological depths via ray tracing (parallelised)
    const size_t nFiltered = bd.voxelIndices.size();
    bd.radDepths.resize(nFiltered, 0.0);

    if (edData) {
        #pragma omp parallel for schedule(dynamic, 256)
        for (size_t i = 0; i < nFiltered; ++i) {
            size_t idx = bd.voxelIndices[i];
            size_t iz = idx / (nx * ny);
            size_t rem = idx % (nx * ny);
            size_t iy = rem / nx;
            size_t ix = rem % nx;
            Vec3 targetLPS = doseGrid.voxelToPatient({
                static_cast<double>(ix), static_cast<double>(iy), static_cast<double>(iz)});

            // Trace from source to targetLPS through CT-grid density
            auto rayDepths = SiddonRayTracer::traceRadDepth(sourcePoint, targetLPS, ctGrid, edData);
            if (!rayDepths.empty()) {
                bd.radDepths[i] = rayDepths.back().second; // cumulative depth at target
            }
        }
    } else {
        // No ED data: assume water (rad depth = geometric depth along beam axis)
        for (size_t i = 0; i < bd.bevCoords.size(); ++i) {
            // In BEV, depth is along + Y axis: depth = SAD + bev_y
            bd.radDepths[i] = std::max(0.0, m_SAD + bd.bevCoords[i][1]);
        }
    }

    // 3. Compute SSDs per ray
    std::vector<Vec3> rayTargets;
    rayTargets.reserve(beam.getNumOfRays());
    for (size_t r = 0; r < beam.getNumOfRays(); ++r) {
        rayTargets.push_back(beam.getRay(r)->getTargetPoint());
    }

    if (edData) {
        bd.ssds = SSDCalculator::computeBeamSSDs(sourcePoint, rayTargets, ctGrid, edData);
    } else {
        bd.ssds.assign(beam.getNumOfRays(), m_SAD); // Default SSD
    }

    // 4. Select nearest SSD kernel and build 2D interpolators
    // Find mean SSD for this beam
    double meanSSD = 0.0;
    int validCount = 0;
    for (double ssd : bd.ssds) {
        if (ssd > 0) { meanSSD += ssd; validCount++; }
    }
    meanSSD = validCount > 0 ? meanSSD / validCount : m_SAD;

    // Find nearest kernel entry by SSD
    size_t bestKernelIdx = 0;
    double bestDist = std::numeric_limits<double>::max();
    for (size_t ki = 0; ki < m_kernelEntries.size(); ++ki) {
        double dist = std::abs(m_kernelEntries[ki].SSD - meanSSD);
        if (dist < bestDist) {
            bestDist = dist;
            bestKernelIdx = ki;
        }
    }

    const auto& selectedKernel = m_kernelEntries[bestKernelIdx];

    // ─── Build 2D kernel interpolators (matching matRad exactly) ───
    // matRad grid convention: -N*res : res : (N-1)*res  → 2N elements (even size)
    double res = m_convResolution;  // 0.5 mm
    double sigma = m_penumbraFWHM / std::sqrt(8.0 * std::log(2.0));

    // Grid limits (in units of res pixels)
    int fieldLimit = static_cast<int>(std::ceil(m_bixelWidth / (2.0 * res)));  // 7
    int gaussLimit = static_cast<int>(std::ceil(5.0 * sigma / res));           // 22
    int kernelLimit = static_cast<int>(std::ceil(m_kernelCutOff / res));       // 359
    int kernelConvLimit = fieldLimit + gaussLimit + kernelLimit;                // 388

    // Sizes (all even, matching matRad's 2*limit convention)
    size_t fieldN = static_cast<size_t>(std::floor(m_bixelWidth / res));       // 14
    size_t gaussConvN = static_cast<size_t>(2 * (fieldLimit + gaussLimit));     // 58
    size_t kernelN = static_cast<size_t>(2 * kernelLimit);                     // 718
    size_t kernelConvN = static_cast<size_t>(2 * kernelConvLimit);             // 776

    // Build Gaussian filter on grid: -gaussLimit*res : res : (gaussLimit-1)*res → 2*gaussLimit elements
    size_t gaussN = static_cast<size_t>(2 * gaussLimit);  // 44
    std::vector<double> gaussFilter(gaussN * gaussN, 0.0);
    double gaussNorm = 1.0 / (2.0 * M_PI * sigma * sigma / (res * res));
    for (size_t ri = 0; ri < gaussN; ++ri) {
        for (size_t ci = 0; ci < gaussN; ++ci) {
            double gx = (static_cast<double>(ri) - gaussLimit) * res;
            double gz = (static_cast<double>(ci) - gaussLimit) * res;
            gaussFilter[ri * gaussN + ci] = gaussNorm * std::exp(-(gx * gx + gz * gz) / (2.0 * sigma * sigma));
        }
    }

    // Build uniform fluence field: ones(fieldN, fieldN)
    std::vector<double> Fpre(fieldN * fieldN, 1.0);

    // Convolve fluence with Gaussian using FFT (pad to gaussConvN×gaussConvN)
    // Matches: Fpre = real(ifft2(fft2(Fpre,gaussConvN,gaussConvN).*fft2(gaussFilter,gaussConvN,gaussConvN)))
    Fpre = FFT2D::convolve2DPadded(Fpre, fieldN, fieldN, gaussFilter, gaussN, gaussN, gaussConvN);

    // Convolution output grid: -kernelConvLimit*res : res : (kernelConvLimit-1)*res
    double convGridMin = -static_cast<double>(kernelConvLimit) * res;
    double convGridMax = (static_cast<double>(kernelConvLimit) - 1.0) * res;

    // For each kernel component, create 2D convolved profile
    for (int k = 0; k < 3; ++k) {
        const auto& kernel1d = (k == 0) ? selectedKernel.kernel1 :
                               (k == 1) ? selectedKernel.kernel2 : selectedKernel.kernel3;

        if (kernel1d.empty() || m_kernelPos.empty()) {
            // Fallback: simple zero interpolator
            std::vector<double> kernData(kernelConvN * kernelConvN, 0.0);
            bd.kernelInterps[k].setGrid(convGridMin, convGridMax, kernelConvN,
                                         convGridMin, convGridMax, kernelConvN, kernData);
            continue;
        }

        // Interpolate 1D radial kernel onto 2D grid: kernelN×kernelN
        // Grid: -kernelLimit*res : res : (kernelLimit-1)*res
        std::vector<double> kernel2d(kernelN * kernelN, 0.0);
        for (size_t ri = 0; ri < kernelN; ++ri) {
            for (size_t ci = 0; ci < kernelN; ++ci) {
                double x = (static_cast<double>(ri) - kernelLimit) * res;
                double z = (static_cast<double>(ci) - kernelLimit) * res;
                double r = std::sqrt(x * x + z * z);

                // Linear interpolation on radial kernel (extrapolate to 0 outside range)
                double kVal = 0.0;
                if (r <= m_kernelPos.back()) {
                    size_t lo = 0;
                    for (size_t pi = 0; pi + 1 < m_kernelPos.size(); ++pi) {
                        if (m_kernelPos[pi + 1] >= r) { lo = pi; break; }
                    }
                    double t = 0.0;
                    if (lo + 1 < m_kernelPos.size() && m_kernelPos[lo + 1] != m_kernelPos[lo]) {
                        t = (r - m_kernelPos[lo]) / (m_kernelPos[lo + 1] - m_kernelPos[lo]);
                    }
                    t = std::clamp(t, 0.0, 1.0);
                    kVal = (1.0 - t) * kernel1d[lo] + t * kernel1d[std::min(lo + 1, kernel1d.size() - 1)];
                }
                kernel2d[ri * kernelN + ci] = kVal;
            }
        }

        // Convolve Fpre with kernel using FFT (pad to kernelConvN×kernelConvN)
        // This matches: convMx = real(ifft2(fft2(Fpre,kernelConvN,kernelConvN).*fft2(kernelMx,kernelConvN,kernelConvN)))
        auto kernelConv = FFT2D::convolve2DPadded(Fpre, gaussConvN, gaussConvN,
                                                    kernel2d, kernelN, kernelN, kernelConvN);

        bd.kernelInterps[k].setGrid(convGridMin, convGridMax, kernelConvN,
                                     convGridMin, convGridMax, kernelConvN, kernelConv);
    }

    return bd;
}

// ────────────────────────────────────────────────
// initRay: select voxels within lateral cutoff and compute iso-plane positions
// ────────────────────────────────────────────────
PencilBeamEngine::RayVoxelData PencilBeamEngine::initRay(
    const Ray& ray, const Beam& beam, const BeamData& beamData)
{
    RayVoxelData rd;

    // Ray target in BEV
    Vec3 rayPosBev = ray.getRayPosBev();

    // For each voxel, check if within lateral cutoff from this ray
    for (size_t i = 0; i < beamData.bevCoords.size(); ++i) {
        const Vec3& bev = beamData.bevCoords[i];
        double geoDist = beamData.geoDistances[i];

        if (geoDist < 1e-6) continue;

        // BEV: x is lateral, y is depth (along beam), z is lateral
        // Lateral offset in BEV from ray position
        double latX = bev[0] - rayPosBev[0];
        double latZ = bev[2] - rayPosBev[2];

        // Project to isocenter plane: lat_iso = lat_bev * SAD / depth
        // Depth in BEV is along y axis: depth = SAD + bev[1] (since source is at -SAD)
        double depth = m_SAD + bev[1];
        if (depth < 1.0) continue; // Behind source

        double isoLatX = latX * m_SAD / depth;
        double isoLatZ = latZ * m_SAD / depth;

        double radialDist = std::sqrt(isoLatX * isoLatX + isoLatZ * isoLatZ);

        if (radialDist <= m_lateralCutOff) {
            rd.localIndices.push_back(i);
            rd.isoLatX.push_back(isoLatX);
            rd.isoLatZ.push_back(isoLatZ);
        }
    }

    return rd;
}

// ────────────────────────────────────────────────
// calcBixelDose: Bortfeld formula for a single bixel
// ────────────────────────────────────────────────
std::vector<double> PencilBeamEngine::calcBixelDose(
    const RayVoxelData& rayData,
    const BeamData& beamData,
    double SAD)
{
    size_t n = rayData.localIndices.size();
    std::vector<double> dose(n, 0.0);

    for (size_t vi = 0; vi < n; ++vi) {
        size_t idx = rayData.localIndices[vi];
        double radDepth = beamData.radDepths[idx];
        double geoDist = beamData.geoDistances[idx];
        double isoX = rayData.isoLatX[vi];
        double isoZ = rayData.isoLatZ[vi];

        if (radDepth <= 0.0 || geoDist < 1.0) continue;

        // Bortfeld depth-dose formula for 3 components
        double totalDose = 0.0;
        for (int k = 0; k < 3; ++k) {
            double beta = m_betas[k];
            double denom = beta - m_mu;

            double depthDose = 0.0;
            if (std::abs(denom) > 1e-12) {
                depthDose = (beta / denom) *
                    (std::exp(-m_mu * radDepth) - std::exp(-beta * radDepth));
            } else {
                // Special case: beta ≈ mu → use L'Hôpital limit = mu * d * exp(-mu*d)
                depthDose = m_mu * radDepth * std::exp(-m_mu * radDepth);
            }

            // Lateral kernel profile (from convolved interpolator)
            double lateralVal = beamData.kernelInterps[k](isoX, isoZ);

            totalDose += depthDose * lateralVal;
        }

        // Inverse square correction
        double invSq = (SAD * SAD) / (geoDist * geoDist);
        totalDose *= invSq;

        // Clamp negative (can arise from FFT ringing)
        if (totalDose < 0.0) totalDose = 0.0;

        dose[vi] = totalDose;
    }

    return dose;
}

// ────────────────────────────────────────────────
// calculateDij: main dij computation loop
// ────────────────────────────────────────────────
DoseInfluenceMatrix PencilBeamEngine::calculateDij(
    const Plan& plan,
    const Stf& stf,
    const PatientData& patientData,
    const Grid& doseGrid)
{
    auto startTime = std::chrono::steady_clock::now();

    m_bixelWidth = plan.getStfProperties().bixelWidth;
    initDoseCalc(plan, doseGrid);

    // Apply OMP thread setting from options
    const auto& opts = m_options;
#ifdef _OPENMP
    if (opts.numThreads > 0) {
        omp_set_num_threads(opts.numThreads);
    }
    Logger::info("PencilBeamEngine: Using " +
        std::to_string(opts.numThreads > 0 ? opts.numThreads : omp_get_max_threads()) +
        " OpenMP threads");
#endif
    Logger::info("PencilBeamEngine: absoluteThreshold=" +
        std::to_string(opts.absoluteThreshold) +
        ", relativeThreshold=" + std::to_string(opts.relativeThreshold));

    // Ensure ED volume exists
    if (!patientData.getEDVolume()) {
        Logger::warn("PencilBeamEngine: No electron density volume. Using unit density.");
    }

    // Collect union of all structure voxel indices (in CT grid), then map to dose grid
    const auto& ctGrid = patientData.getGrid();
    const auto* structSet = patientData.getStructureSet();
    std::unordered_set<size_t> uniqueCtIndices;

    if (structSet) {
        for (size_t s = 0; s < structSet->getCount(); ++s) {
            const auto* structure = structSet->getStructure(s);
            if (!structure) continue;
            const auto& voxels = structure->getVoxelIndices();
            uniqueCtIndices.insert(voxels.begin(), voxels.end());
        }
    }

    std::vector<size_t> ctIndicesVec(uniqueCtIndices.begin(), uniqueCtIndices.end());

    // Map to dose grid voxel indices
    std::vector<size_t> doseVoxelIndices;
    if (!ctIndicesVec.empty()) {
        doseVoxelIndices = Grid::resampleMaskNearestToGrid(ctGrid, doseGrid, ctIndicesVec);
    }

    if (doseVoxelIndices.empty()) {
        // Fallback: use all voxels in dose grid
        Logger::warn("PencilBeamEngine: No structure voxels found. Using all dose grid voxels.");
        doseVoxelIndices.resize(doseGrid.getNumVoxels());
        std::iota(doseVoxelIndices.begin(), doseVoxelIndices.end(), 0);
    }

    // Count total bixels across all beams
    size_t totalBixels = stf.getTotalNumOfBixels();
    size_t numDoseVoxels = doseGrid.getNumVoxels();

    Logger::info("PencilBeamEngine: Dij dimensions = " +
        std::to_string(numDoseVoxels) + " voxels x " +
        std::to_string(totalBixels) + " bixels");
    Logger::info("PencilBeamEngine: Active voxels (in structures) = " +
        std::to_string(doseVoxelIndices.size()));

    // Allocate Dij
    DoseInfluenceMatrix dij(numDoseVoxels, totalBixels);

    // Precompute bixel offsets for each beam and ray
    // bixelOffsets[bi] = global bixel index where beam bi starts
    size_t numBeams = stf.getCount();
    std::vector<size_t> beamBixelOffset(numBeams + 1, 0);
    for (size_t bi = 0; bi < numBeams; ++bi) {
        const auto* beam = stf.getBeam(bi);
        size_t beamBixels = 0;
        if (beam) {
            for (size_t ri = 0; ri < beam->getNumOfRays(); ++ri) {
                const auto* ray = beam->getRay(ri);
                if (ray) beamBixels += ray->getNumOfBixels();
            }
        }
        beamBixelOffset[bi + 1] = beamBixelOffset[bi] + beamBixels;
    }

    // Stats
    size_t totalEntriesBefore = 0;
    size_t totalEntriesAfter = 0;

    // Main loop: beams → parallel rays → bixels
    for (size_t bi = 0; bi < numBeams; ++bi) {
        if (m_cancelFlag && m_cancelFlag->load()) {
            Logger::warn("PencilBeamEngine: Cancelled by user.");
            break;
        }

        const auto* beam = stf.getBeam(bi);
        if (!beam) continue;

        size_t numRays = beam->getNumOfRays();

        Logger::info("PencilBeamEngine: Processing beam " + std::to_string(bi + 1) +
            "/" + std::to_string(numBeams) +
            " (gantry=" + std::to_string(beam->getGantryAngle()) +
            " deg, couch=" + std::to_string(beam->getCouchAngle()) +
            " deg, " + std::to_string(numRays) + " rays)");

        if (m_progressCallback) {
            m_progressCallback(static_cast<int>(bi), static_cast<int>(numBeams),
                "Beam " + std::to_string(bi + 1) + "/" + std::to_string(numBeams));
        }

        // Per-beam initialization (parallelised radDepth + BEV)
        BeamData beamData = initBeam(*beam, patientData, doseGrid, doseVoxelIndices);

        // Precompute per-ray bixel offset within this beam
        std::vector<size_t> rayBixelOff(numRays + 1, 0);
        for (size_t ri = 0; ri < numRays; ++ri) {
            const auto* ray = beam->getRay(ri);
            rayBixelOff[ri + 1] = rayBixelOff[ri] + (ray ? ray->getNumOfBixels() : 0);
        }

        size_t beamStart = beamBixelOffset[bi];
        size_t beamEntriesBefore = 0;
        size_t beamEntriesAfter = 0;

        // ── Parallel ray processing ──
        // Each thread collects COO entries locally, then merges under critical section.
        #pragma omp parallel reduction(+:beamEntriesBefore,beamEntriesAfter)
        {
            // Thread-local COO buffers
            std::vector<size_t> localRows;
            std::vector<size_t> localCols;
            std::vector<double> localVals;

            #pragma omp for schedule(dynamic, 4)
            for (size_t ri = 0; ri < numRays; ++ri) {
                const auto* ray = beam->getRay(ri);
                if (!ray) continue;

                // Per-ray initialization (voxel selection within lateral cutoff)
                RayVoxelData rayData = initRay(*ray, *beam, beamData);

                if (rayData.localIndices.empty()) continue;

                // Per-bixel dose calculation (photons: 1 bixel per ray)
                auto bixelDose = calcBixelDose(rayData, beamData, m_SAD);

                size_t bixelCol = beamStart + rayBixelOff[ri];

                // Find max dose for this bixel (for relative threshold)
                double maxDose = 0.0;
                for (size_t vi = 0; vi < bixelDose.size(); ++vi) {
                    if (bixelDose[vi] > maxDose) maxDose = bixelDose[vi];
                }

                double relThresh = maxDose * opts.relativeThreshold;
                double threshold = std::max(relThresh, opts.absoluteThreshold);

                beamEntriesBefore += bixelDose.size();

                // Filter and collect entries
                for (size_t vi = 0; vi < rayData.localIndices.size(); ++vi) {
                    if (bixelDose[vi] >= threshold) {
                        size_t voxelIdx = beamData.voxelIndices[rayData.localIndices[vi]];
                        localRows.push_back(voxelIdx);
                        localCols.push_back(bixelCol);
                        localVals.push_back(bixelDose[vi]);
                        beamEntriesAfter++;
                    }
                }
            }

            // Merge thread-local COO into shared Dij
            #pragma omp critical
            {
                dij.appendBatch(localRows, localCols, localVals);
            }
        } // end parallel

        totalEntriesBefore += beamEntriesBefore;
        totalEntriesAfter += beamEntriesAfter;
        
        double pctKept = beamEntriesBefore > 0
            ? 100.0 * static_cast<double>(beamEntriesAfter) / static_cast<double>(beamEntriesBefore)
            : 0.0;
        Logger::info("PencilBeamEngine: Beam " + std::to_string(bi + 1) +
            ": " + std::to_string(beamEntriesAfter) + " entries kept" +
            " (threshold removed " + std::to_string(100.0 - pctKept).substr(0, 5) + "%" +
            ", COO size ~" + std::to_string(dij.getNumNonZeros() * 24 / (1024 * 1024)) + " MB)");
    }

    // Finalize: convert to CSR sparse format
    Logger::info("PencilBeamEngine: Finalizing Dij (sparse conversion)...");
    Logger::info("PencilBeamEngine: Total entries before threshold: " +
        std::to_string(totalEntriesBefore) +
        ", after: " + std::to_string(totalEntriesAfter) +
        " (" + std::to_string(totalEntriesBefore > 0
            ? 100.0 * static_cast<double>(totalEntriesAfter) / static_cast<double>(totalEntriesBefore)
            : 0.0).substr(0, 5) + "% kept)");
    dij.finalize();

    auto endTime = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(endTime - startTime).count();
    Logger::info("PencilBeamEngine: Dij computation complete in " +
        std::to_string(elapsed) + "s" +
        " (nnz=" + std::to_string(dij.getNumNonZeros()) + ")");

    if (m_progressCallback) {
        m_progressCallback(static_cast<int>(numBeams), static_cast<int>(numBeams), "Done");
    }

    return dij;
}

// ────────────────────────────────────────────────
// calculateDose: forward dose from Dij and weights
// ────────────────────────────────────────────────
DoseMatrix PencilBeamEngine::calculateDose(
    const DoseInfluenceMatrix& dij,
    const std::vector<double>& weights,
    const Grid& grid)
{
    DoseMatrix dose;
    dose.setGrid(grid);
    dose.allocate();

    auto doseVec = dij.computeDose(weights);

    // Copy into dose matrix
    auto dims = grid.getDimensions();
    size_t nx = dims[0], ny = dims[1], nz = dims[2];

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                size_t flatIdx = i + j * nx + k * nx * ny;
                if (flatIdx < doseVec.size()) {
                    dose.at(i, j, k) = doseVec[flatIdx];
                }
            }
        }
    }

    return dose;
}

} // namespace optirad
