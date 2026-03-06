#pragma once

#include "../IDoseEngine.hpp"
#include "../GridInterpolant2D.hpp"
#include "../FFT2D.hpp"
#include "core/Machine.hpp"
#include <vector>
#include <array>

namespace optirad {

/**
 * Photon Pencil Beam SVD Dose Engine.
 *
 * Implements the Bortfeld (1993) collapsed-cone / SVD pencil beam algorithm
 * as used in matRad's matRad_PhotonPencilBeamSVDEngine.
 *
 * The dose from a single bixel at radiological depth d_rad is:
 *   D_k(d) = beta_k / (beta_k - mu) * (exp(-mu * d_rad) - exp(-beta_k * d_rad))
 *   D_total = sum_k(D_k * K_k(x_iso, z_iso)) * (SAD / d_geo)^2
 *
 * Where:
 *   - beta_k: attenuation coefficients for 3 kernel components (primary, 1st scatter, 2nd scatter)
 *   - mu (m): primary absorption coefficient
 *   - K_k: lateral kernel profiles pre-convolved with fluence and gaussian via FFT
 *   - SAD: source-axis distance
 *   - d_geo: geometric distance from source to voxel
 */
class PencilBeamEngine : public IDoseEngine {
public:
    std::string getName() const override;

    DoseInfluenceMatrix calculateDij(
        const Plan& plan,
        const Stf& stf,
        const PatientData& patientData,
        const Grid& doseGrid) override;

    DoseMatrix calculateDose(
        const DoseInfluenceMatrix& dij,
        const std::vector<double>& weights,
        const Grid& grid) override;

private:
    // ── Initialization ──
    void initDoseCalc(const Plan& plan, const Grid& doseGrid);

    // ── Per-beam setup ──
    struct BeamData {
        std::vector<size_t> voxelIndices;    // Dose grid voxel indices affected by this beam
        std::vector<double> radDepths;       // Radiological depths at each voxel
        std::vector<double> geoDistances;    // Geometric distances from source
        std::vector<Vec3>   bevCoords;       // BEV coordinates of each voxel
        std::vector<double> ssds;            // SSD per ray
        GridInterpolant2D kernelInterps[3];  // Convolved kernel interpolators (3 components)
    };

    BeamData initBeam(const Beam& beam, const PatientData& patientData,
                      const Grid& doseGrid, const std::vector<size_t>& allVoxelIndices);

    // ── Per-ray: compute lateral distances and voxel selection ──
    struct RayVoxelData {
        std::vector<size_t> localIndices;    // Indices into BeamData arrays for voxels near this ray
        std::vector<double> isoLatX;         // Lateral X distance at isocenter plane
        std::vector<double> isoLatZ;         // Lateral Z distance at isocenter plane
    };

    RayVoxelData initRay(const Ray& ray, const Beam& beam, const BeamData& beamData);

    // ── Per-bixel dose calculation (Bortfeld formula) ──
    std::vector<double> calcBixelDose(
        const RayVoxelData& rayData,
        const BeamData& beamData,
        double SAD);

    // ── Machine parameters (cached from init) ──
    std::array<double, 3> m_betas;
    double m_mu;  // primary absorption coefficient
    double m_SAD;
    double m_SCD;
    double m_penumbraFWHM;
    double m_geometricLateralCutOff = 50.0; // mm, geometric cutoff for voxel filtering (matRad default)
    double m_kernelCutOff = 179.5; // mm, full kernel extent (for convolution grid)
    double m_lateralCutOff = 50.0; // mm, effective lateral cutoff used in ray voxel filtering
    double m_bixelWidth;

    // Kernel data from machine
    std::vector<double> m_kernelPos; // radial positions
    std::vector<KernelEntry> m_kernelEntries; // kernels at various SSDs

    // Convolution parameters
    double m_convResolution = 0.5; // mm, internal convolution grid resolution
};

} // namespace optirad
