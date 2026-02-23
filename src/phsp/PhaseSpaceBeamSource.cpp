#include "phsp/PhaseSpaceBeamSource.hpp"
#include "phsp/IAEAPhspReader.hpp"
#include "utils/Logger.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace optirad {

void PhaseSpaceBeamSource::configure(const Machine& machine,
                                      double gantryAngle,
                                      double collimatorAngle,
                                      double couchAngle,
                                      const std::array<double, 3>& isocenter) {
    m_machine = machine;
    m_gantryAngle = gantryAngle;
    m_collimatorAngle = collimatorAngle;
    m_couchAngle = couchAngle;
    m_isocenter = isocenter;

    // Compute rotation matrix (gantry + couch)
    m_rotMatrix = getRotationMatrix(gantryAngle, couchAngle);

    // Compute source position
    computeSourcePosition();

    m_built = false;
}

void PhaseSpaceBeamSource::computeSourcePosition() {
    double SAD = m_machine.getSAD(); // mm
    // Source in BEV is at (0, -SAD, 0); transform to LPS
    Vec3 sourceBEV = {0.0, -SAD, 0.0};
    Vec3 sourceLPS = vecMulMatTranspose(sourceBEV, m_rotMatrix);
    m_sourcePos = vecAdd(sourceLPS, m_isocenter);
}

void PhaseSpaceBeamSource::build(int64_t maxParticles, int64_t vizSampleSize) {
    const auto& geom = m_machine.getGeometry();

    if (geom.phaseSpaceFileNames.empty()) {
        throw std::runtime_error("PhaseSpaceBeamSource: no phase-space files configured");
    }

    // Parse header from first file (structure assumed identical across all files)
    std::string firstHeaderPath = geom.phaseSpaceDir + "/" +
                                   geom.phaseSpaceFileNames[0] + ".IAEAheader";
    m_headerInfo = IAEAHeaderParser::parse(firstHeaderPath);

    const int numFiles = static_cast<int>(geom.phaseSpaceFileNames.size());
    Logger::info("PhaseSpaceBeamSource: building from " +
                 std::to_string(numFiles) + " PSF file(s), " +
                 "gantry=" + std::to_string(m_gantryAngle) + "°");

    // Determine total particles to read, spread across all files
    int64_t readCount = maxParticles;
    if (readCount <= 0) {
        // Default: up to 5M total, drawn evenly from all files
        readCount = std::min(m_headerInfo.totalParticles * static_cast<int64_t>(numFiles),
                             int64_t(5000000));
    }

    // Distribute particles evenly across PSF files for balanced MC statistics
    int64_t perFileCount = readCount / numFiles;
    int64_t remainder    = readCount % numFiles;

    // Read all PSF files in parallel (OpenMP)
    std::vector<PhaseSpaceData> fileData(numFiles);

    #pragma omp parallel for schedule(static)
    for (int f = 0; f < numFiles; ++f) {
        int64_t thisCount = perFileCount + (f < static_cast<int>(remainder) ? 1 : 0);
        std::string phspPath = geom.phaseSpaceDir + "/" +
                               geom.phaseSpaceFileNames[f] + ".IAEAphsp";

        if (thisCount < m_headerInfo.totalParticles / 2) {
            fileData[f] = IAEAPhspReader::readSampled(phspPath, m_headerInfo, thisCount);
        } else {
            fileData[f] = IAEAPhspReader::readSubset(phspPath, m_headerInfo, 0, thisCount);
        }
    }

    // Merge results from all files
    size_t totalParticles = 0;
    for (int f = 0; f < numFiles; ++f) {
        totalParticles += fileData[f].size();
    }
    m_data.particles().reserve(totalParticles);
    for (int f = 0; f < numFiles; ++f) {
        auto& src = fileData[f].particles();
        m_data.particles().insert(m_data.particles().end(),
                                   std::make_move_iterator(src.begin()),
                                   std::make_move_iterator(src.end()));
    }

    Logger::debug("PhaseSpaceBeamSource: read " + std::to_string(m_data.size()) + " particles");

    // Apply aperture filtering (jaws)
    double jawX1 = geom.defaultFieldSize[0] / 2.0; // half-width in mm
    double jawX2 = geom.defaultFieldSize[0] / 2.0;
    double jawY1 = geom.defaultFieldSize[1] / 2.0;
    double jawY2 = geom.defaultFieldSize[1] / 2.0;
    double scoringPlaneZ = (m_headerInfo.zRange[0] + m_headerInfo.zRange[1]) / 2.0; // cm

    size_t beforeFilter = m_data.size();
    m_data.filterByJaws(jawX1, jawX2, jawY1, jawY2,
                        m_machine.getSAD(), scoringPlaneZ);
    Logger::debug("PhaseSpaceBeamSource: after jaw filter: " + std::to_string(m_data.size()) +
                  " / " + std::to_string(beforeFilter) + " particles");

    // Transform to patient coordinates
    transformToPatientCoords();

    // Compute metrics
    m_metrics = m_data.computeMetrics();

    // Create visualization sample
    if (vizSampleSize > 0 && static_cast<int64_t>(m_data.size()) > vizSampleSize) {
        m_vizSample = m_data.sample(static_cast<size_t>(vizSampleSize));
    } else {
        m_vizSample = m_data;
    }

    m_built = true;

    Logger::debug("PhaseSpaceBeamSource: built successfully, " +
                 std::to_string(m_metrics.totalCount) + " particles, " +
                 "mean E=" + std::to_string(m_metrics.meanEnergy) + " MeV, " +
                 "viz sample=" + std::to_string(m_vizSample.size()));
}

void PhaseSpaceBeamSource::transformToPatientCoords() {
    // Transform from PSF reference frame to patient LPS coordinates.
    //
    // PSF frame (from IAEA header):
    //   Origin: target bottom center
    //   +X: toward linac (parallel to upper jaws)
    //   +Z: toward patient
    //
    // Our BEV convention:
    //   Source at (0, -SAD, 0)
    //   Isocenter at (0, 0, 0)
    //   Rays project from source through isocenter toward +Y
    //
    // Mapping PSF → BEV:
    //   PSF_X → BEV_X (lateral)
    //   PSF_Y → BEV_Z (lateral, perpendicular to jaws)
    //   PSF_Z → BEV_Y + offset (along beam axis, toward patient)
    //   The PSF scoring plane is at some Z above the isocenter;
    //   we need to account for the scoring plane distance from source.
    //
    // Then BEV → LPS via rotation matrix + isocenter offset.

    double SAD = m_machine.getSAD(); // mm

    for (auto& p : m_data.particles()) {
        // Current position is in mm (converted during read), directions are cosines
        // PSF coords: (px, py, pz) where pz is along beam toward patient

        // Map PSF → BEV:
        //   BEV_x = PSF_x
        //   BEV_z = PSF_y  (the second lateral dimension)
        //   BEV_y = PSF_z - SAD  (shift so isocenter is at y=0)
        // PSF_z is the scoring plane Z in mm (already converted from cm)
        Vec3 posBEV = {p.position[0], p.position[2] - SAD, p.position[1]};

        // Direction cosines: similarly remap
        Vec3 dirBEV = {p.direction[0], p.direction[2], p.direction[1]};

        // Apply collimator rotation (rotation around BEV Y-axis, i.e. beam axis)
        if (std::abs(m_collimatorAngle) > 0.01) {
            double rad = m_collimatorAngle * M_PI / 180.0;
            double cc = std::cos(rad), sc = std::sin(rad);
            double newX = cc * posBEV[0] - sc * posBEV[2];
            double newZ = sc * posBEV[0] + cc * posBEV[2];
            posBEV[0] = newX;
            posBEV[2] = newZ;

            double newDirX = cc * dirBEV[0] - sc * dirBEV[2];
            double newDirZ = sc * dirBEV[0] + cc * dirBEV[2];
            dirBEV[0] = newDirX;
            dirBEV[2] = newDirZ;
        }

        // BEV → LPS via gantry+couch rotation matrix
        Vec3 posLPS = vecMulMatTranspose(posBEV, m_rotMatrix);
        Vec3 dirLPS = vecMulMatTranspose(dirBEV, m_rotMatrix);

        // Translate to isocenter position
        p.position = vecAdd(posLPS, m_isocenter);
        p.direction = dirLPS;
    }

    // Do the same for viz sample if it's been populated
    // (not yet at this point — it's populated after this call)
}

std::vector<std::pair<double, int64_t>> PhaseSpaceBeamSource::computeEnergyHistogram(int numBins) const {
    std::vector<std::pair<double, int64_t>> histogram;
    if (m_data.empty() || numBins <= 0) return histogram;

    double minE = m_metrics.minEnergy;
    double maxE = m_metrics.maxEnergy;
    if (maxE <= minE) {
        histogram.push_back({minE, static_cast<int64_t>(m_data.size())});
        return histogram;
    }

    double binWidth = (maxE - minE) / numBins;
    histogram.resize(numBins, {0.0, 0});

    for (int i = 0; i < numBins; ++i) {
        histogram[i].first = minE + (i + 0.5) * binWidth; // bin center
    }

    for (const auto& p : m_data.particles()) {
        int bin = static_cast<int>((p.energy - minE) / binWidth);
        bin = std::clamp(bin, 0, numBins - 1);
        histogram[bin].second++;
    }

    return histogram;
}

} // namespace optirad
