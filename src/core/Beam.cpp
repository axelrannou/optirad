#include "Beam.hpp"
#include <cmath>
#include <algorithm>

namespace optirad {

void Beam::setGantryAngle(double angle) { m_gantryAngle = angle; }
void Beam::setCouchAngle(double angle) { m_couchAngle = angle; }

void Beam::setIsocenter(double x, double y, double z) {
    m_isocenter = {x, y, z};
}

void Beam::setIsocenter(const Vec3& iso) {
    m_isocenter = iso;
}

double Beam::getGantryAngle() const { return m_gantryAngle; }
double Beam::getCouchAngle() const { return m_couchAngle; }

void Beam::addRay(const Ray& ray) {
    m_rays.push_back(ray);
}

void Beam::addRay(Ray&& ray) {
    m_rays.push_back(std::move(ray));
}

const Ray* Beam::getRay(size_t index) const {
    if (index >= m_rays.size()) return nullptr;
    return &m_rays[index];
}

Ray* Beam::getRay(size_t index) {
    if (index >= m_rays.size()) return nullptr;
    return &m_rays[index];
}

std::vector<size_t> Beam::getNumOfBixelsPerRay() const {
    std::vector<size_t> counts;
    counts.reserve(m_rays.size());
    for (const auto& ray : m_rays) {
        counts.push_back(ray.getNumOfBixels());
    }
    return counts;
}

size_t Beam::getTotalNumOfBixels() const {
    size_t total = 0;
    for (const auto& ray : m_rays) {
        total += ray.getNumOfBixels();
    }
    return total;
}

void Beam::computeSourcePoints() {
    // Source point in BEV: source is at (0, -SAD, 0)
    m_sourcePointBev = {0.0, -m_SAD, 0.0};

    // Get the rotation matrix and its transpose for row-vector multiplication
    Mat3 rotMat = getRotationMatrix(m_gantryAngle, m_couchAngle);
    Mat3 rotMatT = transpose(rotMat);

    // Source point in LPS: sourcePoint_bev * rotMat^T
    m_sourcePoint = vecMulMatTranspose(m_sourcePointBev, rotMat);
}

void Beam::initRaysFromPositions(const std::vector<Vec3>& rayPositionsBev) {
    m_rays.clear();
    m_rays.reserve(rayPositionsBev.size());

    // Compute source point first
    computeSourcePoints();

    // Get rotation matrix for BEV->LPS transformation
    Mat3 rotMat = getRotationMatrix(m_gantryAngle, m_couchAngle);

    for (const auto& posBev : rayPositionsBev) {
        Ray ray;

        // Ray position in BEV
        ray.setRayPosBev(posBev);

        // Target point in BEV: target = (2*x, SAD, 2*z) 
        // This projects the ray through the isocenter to a point at distance SAD on the other side
        Vec3 targetBev = {2.0 * posBev[0], m_SAD, 2.0 * posBev[2]};
        ray.setTargetPointBev(targetBev);

        // Transform to LPS using rotation matrix transpose (row-vector convention)
        Vec3 posLps = vecMulMatTranspose(posBev, rotMat);
        Vec3 targetLps = vecMulMatTranspose(targetBev, rotMat);

        ray.setRayPos(posLps);
        ray.setTargetPoint(targetLps);

        m_rays.push_back(std::move(ray));
    }
}

void Beam::computePhotonRayCorners() {
    Mat3 rotMat = getRotationMatrix(m_gantryAngle, m_couchAngle);
    double halfBW = m_bixelWidth / 2.0;
    double scdRatio = m_SCD / m_SAD;

    for (auto& ray : m_rays) {
        const Vec3& pos = ray.getRayPosBev();

        // 4 corners of the beamlet at the isocenter plane (in BEV, then rotated to LPS)
        // Corner order: (+x,+z), (-x,+z), (-x,-z), (+x,-z)
        std::array<Vec3, 4> cornersBev = {{
            {pos[0] + halfBW, 0.0, pos[2] + halfBW},
            {pos[0] - halfBW, 0.0, pos[2] + halfBW},
            {pos[0] - halfBW, 0.0, pos[2] - halfBW},
            {pos[0] + halfBW, 0.0, pos[2] - halfBW}
        }};

        // Transform corners to LPS
        std::array<Vec3, 4> cornersIso;
        for (int c = 0; c < 4; ++c) {
            cornersIso[c] = vecMulMatTranspose(cornersBev[c], rotMat);
        }
        ray.setBeamletCornersAtIso(cornersIso);

        // Corners at SCD plane: offset = (0, SCD-SAD, 0), scaled position = SCD/SAD * corner_bev
        // rayCorners_SCD = (repmat([0, SCD-SAD, 0], 4, 1) + (SCD/SAD) * corner_bev) * rotMat^T
        Vec3 scdOffset = {0.0, m_SCD - m_SAD, 0.0};
        std::array<Vec3, 4> cornersSCD;
        for (int c = 0; c < 4; ++c) {
            Vec3 scaledCorner = vecScale(cornersBev[c], scdRatio);
            Vec3 scdBev = vecAdd(scdOffset, scaledCorner);
            cornersSCD[c] = vecMulMatTranspose(scdBev, rotMat);
        }
        ray.setRayCornersSCD(cornersSCD);
    }
}

void Beam::setAllRayEnergies(double energy) {
    for (auto& ray : m_rays) {
        ray.setEnergy(energy);
    }
}

void Beam::generateRays(double bixelWidth, const std::array<double, 2>& fieldSize) {
    m_bixelWidth = bixelWidth;
    m_fieldSize = fieldSize;

    // Generate ray positions on a regular grid centered at isocenter (BEV)
    double halfW = fieldSize[0] / 2.0;
    double halfH = fieldSize[1] / 2.0;

    std::vector<Vec3> positions;
    for (double z = -halfH; z <= halfH + 1e-10; z += bixelWidth) {
        for (double x = -halfW; x <= halfW + 1e-10; x += bixelWidth) {
            // Snap to bixel grid
            double snappedX = bixelWidth * std::round(x / bixelWidth);
            double snappedZ = bixelWidth * std::round(z / bixelWidth);
            positions.push_back({snappedX, 0.0, snappedZ});
        }
    }

    initRaysFromPositions(positions);
    computePhotonRayCorners();
    setAllRayEnergies(6.0); // Default photon energy
}

void Beam::generateRaysFromTarget(const std::vector<Vec3>& targetWorldCoords,
                                  double bixelWidth,
                                  const Vec3& ctResolution) {
    m_bixelWidth = bixelWidth;

    // Get the rotation matrix for BEV transformation
    // R performs BEV→LPS (active rotation). For LPS→BEV we need R^T.
    Mat3 rotMat = getRotationMatrix(m_gantryAngle, m_couchAngle);
    Mat3 rotMatInv = transpose(rotMat);

    // Step 1: Subtract isocenter from target world coordinates, rotate into BEV,
    //         and project onto isocenter plane with SAD divergence correction
    std::vector<Vec3> projectedCoords(targetWorldCoords.size());

#ifdef _OPENMP
    #pragma omp parallel for if(targetWorldCoords.size() > 1000)
#endif
    for (size_t idx = 0; idx < targetWorldCoords.size(); ++idx) {
        const auto& worldPos = targetWorldCoords[idx];
        // Correct for isocenter position: isocenter becomes (0,0,0)
        Vec3 isoCoords = vecSub(worldPos, m_isocenter);

        // Rotate from LPS into BEV: v_bev = R^T * v_lps
        Vec3 bev = rotMatInv * isoCoords;

        // SAD divergence correction: project onto isocenter plane
        double scale = m_SAD / (m_SAD + bev[1]);
        Vec3 projected = {bev[0] * scale, 0.0, bev[2] * scale};

        projectedCoords[idx] = projected;
    }

    // Step 2: Snap to bixel grid
    std::vector<Vec3> snappedPositions(projectedCoords.size());
    
#ifdef _OPENMP
    #pragma omp parallel for if(projectedCoords.size() > 1000)
#endif
    for (size_t idx = 0; idx < projectedCoords.size(); ++idx) {
        const auto& p = projectedCoords[idx];
        snappedPositions[idx] = {
            bixelWidth * std::round(p[0] / bixelWidth),
            0.0,
            bixelWidth * std::round(p[2] / bixelWidth)
        };
    }

    // Step 3: Remove duplicates first, then pad if CT resolution is coarser than bixel width
    auto vecLess = [](const Vec3& a, const Vec3& b) {
        if (std::abs(a[0] - b[0]) > 1e-10) return a[0] < b[0];
        if (std::abs(a[2] - b[2]) > 1e-10) return a[2] < b[2];
        return false;
    };
    auto vecEqual = [](const Vec3& a, const Vec3& b) {
        return std::abs(a[0] - b[0]) < 1e-10 &&
               std::abs(a[2] - b[2]) < 1e-10;
    };
    std::sort(snappedPositions.begin(), snappedPositions.end(), vecLess);
    snappedPositions.erase(
        std::unique(snappedPositions.begin(), snappedPositions.end(), vecEqual),
        snappedPositions.end()
    );

    // Step 4: Pad ray positions if CT resolution is coarser than bixel width
    double maxCtRes = std::max({ctResolution[0], ctResolution[1], ctResolution[2]});
    if (bixelWidth < maxCtRes) {
        std::vector<Vec3> padded = snappedPositions;
        int pad = static_cast<int>(std::floor(maxCtRes / bixelWidth));
        for (const auto& pos : snappedPositions) {
            for (int j = -pad; j <= pad; ++j) {
                for (int k = -pad; k <= pad; ++k) {
                    if (j == 0 && k == 0) continue;
                    padded.push_back({
                        pos[0] + j * bixelWidth,
                        0.0,
                        pos[2] + k * bixelWidth
                    });
                }
            }
        }
        snappedPositions = std::move(padded);

        // Remove duplicates again after padding
        std::sort(snappedPositions.begin(), snappedPositions.end(), vecLess);
        snappedPositions.erase(
            std::unique(snappedPositions.begin(), snappedPositions.end(), vecEqual),
            snappedPositions.end()
        );
    }

    // Step 5: Create rays from the unique BEV positions
    initRaysFromPositions(snappedPositions);
    computePhotonRayCorners();
}

} // namespace optirad
