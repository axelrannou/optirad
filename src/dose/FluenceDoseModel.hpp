#pragma once

#include "core/Beam.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace optirad {

class CudaFluenceBackend;

enum class DoseComputeBackend { CPU, Auto, CUDA };

struct FluenceModelOptions {
    /// Fluence pixel size in the isocenter plane (mm); snapped to an odd multiple of the 0.5 mm kernel grid.
    double pixelMm = 1.5;
    DoseComputeBackend backend = DoseComputeBackend::Auto;
    size_t gpuMemoryBudgetMiB = 0;
};

/// Matrix-free pencil beam dose: no Dij is stored.
///
/// The Bortfeld engine is separable, d(v) = sum_k D_k(radDepth_v) * (F (*) K_k)(x_v, z_v), where F is a beam's
/// fluence image and K_k the lateral kernel. Per beam only the voxel geometry is kept (24 B per voxel); the
/// convolution is done with an FFT on every evaluation. This requires divergent bixels
/// (DoseCalcOptions::divergentBixels), the only geometry in which the lateral kernel is shift-invariant.
///
/// Fluence images are n x n floats, row = beam-space Z, column = beam-space X, pixel (n/2, n/2) on the beam axis.
/// Built by PencilBeamEngine::buildFluenceModel(); const methods are thread-safe.
class FluenceDoseModel {
public:
    struct VoxelRecord {
        uint32_t voxel;       // dose-grid voxel index
        float fx, fz;         // pixel coordinates of the voxel projected onto the isocenter plane
        float d[3];           // depth-dose of each kernel component times the inverse-square factor
    };

    FluenceDoseModel(size_t numVoxels, size_t imageSize, double pixelMm, double bixelBoxMm);

    /// Lateral kernel tables on the pixel lattice, (tableSize x tableSize), tableSize odd, centred.
    void addKernel(size_t kernelIdx, const std::array<std::vector<double>, 3>& table, size_t tableSize);
    bool hasKernel(size_t kernelIdx) const { return m_kernelSlot.count(kernelIdx) != 0; }
    void addBeam(size_t kernelIdx, std::vector<VoxelRecord> records);

    size_t numBeams() const { return m_beams.size(); }
    size_t numVoxels() const { return m_numVoxels; }
    size_t imageSize() const { return m_n; }
    size_t imagePixels() const { return m_n * m_n; }
    double pixelMm() const { return m_pixel; }
    size_t memoryBytes() const;
    bool enableCuda(size_t memoryBudgetBytes = 0, bool required = false);
    bool usesCuda() const { return m_cuda != nullptr; }
    size_t cudaMemoryBytes() const;
    bool cudaGeometryResident() const;

    /// Pixel coordinate (fractional) of a BEV isocenter-plane position in mm.
    double toPixel(double mm) const { return mm / m_pixel + static_cast<double>(m_n / 2); }

    /// Fluence image of a beam's bixel weights: every bixel is a box of the bixel width, weighted by its
    /// area coverage of each pixel. `weights` points at the beam's first bixel.
    void rasterizeBixels(const Beam& beam, const double* weights, std::vector<float>& image) const;

    /// Adjoint of rasterizeBixels: out[bixel] += sum over pixels of area coverage x image. `out` points at the
    /// beam's first bixel.
    void gatherBixels(const Beam& beam, const std::vector<float>& image, double* out) const;

    /// Voxels (rows) the model covers; empty = all.
    void setComputedRows(std::vector<uint8_t> mask) { m_rowMask = std::move(mask); }
    const std::vector<uint8_t>& computedRows() const { return m_rowMask; }

    /// dose[v] += sum of the beam's dose at voxel v. `fluence` has imagePixels() entries.
    void addBeamDose(size_t beam, const std::vector<float>& fluence, double* dose) const;

    /// Adjoint of addBeamDose: fluenceGrad = A^T voxelGrad (imagePixels() entries).
    void beamGradient(size_t beam, const double* voxelGrad, std::vector<float>& fluenceGrad) const;

    /// Dose of all beams (OpenMP over beams); an empty fluence image skips the beam.
    std::vector<double> computeDose(const std::vector<std::vector<float>>& fluence) const;

    /// Gradient images of all beams for a voxel gradient (OpenMP over beams).
    std::vector<std::vector<float>> computeGradient(const std::vector<double>& voxelGrad) const;

private:
    friend class CudaFluenceBackend;
    friend class FluenceDoseOperator;
    std::shared_ptr<CudaFluenceBackend> m_cuda;
    struct Spectrum { std::array<std::vector<double>, 3> re, im; };

    size_t m_numVoxels;
    size_t m_n;
    double m_pixel;
    double m_boxMm;
    std::unordered_map<size_t, size_t> m_kernelSlot;
    std::vector<Spectrum> m_spectra;
    struct BeamRecords { size_t kernelSlot; std::vector<VoxelRecord> records; };
    std::vector<BeamRecords> m_beams;
    std::vector<uint8_t> m_rowMask;
};

} // namespace optirad
