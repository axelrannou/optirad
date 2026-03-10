# Dose Calculation Module (`optirad_dose`)

The dose module provides the dose calculation engine, sparse dose influence matrix, ray tracing algorithms, plan analysis, and supporting mathematical tools.

**Library:** `optirad_dose`  
**Dependencies:** `optirad_core`, `optirad_geometry`  
**Optional:** OpenMP (parallel dose calculation)

## Files

| File | Description |
|------|-------------|
| `IDoseEngine.hpp` | Abstract dose engine interface |
| `DoseEngineFactory.hpp/cpp` | Factory for creating engines |
| `engines/PencilBeamEngine.hpp/cpp` | Bortfeld SVD pencil beam engine |
| `DoseInfluenceMatrix.hpp/cpp` | Sparse Dij matrix (COO → CSR) |
| `DoseMatrix.hpp/cpp` | 3D dose result grid |
| `DijSerializer.hpp/cpp` | Binary serialization + caching |
| `SiddonRayTracer.hpp/cpp` | Siddon 1985 ray-voxel tracing |
| `SSDCalculator.hpp/cpp` | Source-to-surface distance |
| `RadDepthCalculator.hpp/cpp` | Radiological depth calculation |
| `PlanAnalysis.hpp/cpp` | DVH, statistics, quality indices |
| `DoseCalcOptions.hpp` | Dose calculation parameters |
| `FFT2D.hpp` | 2D FFT convolution (header-only) |
| `GridInterpolant2D.hpp` | 2D bilinear interpolation (header-only) |

## Interfaces

### IDoseEngine

Abstract interface for dose calculation engines.

```cpp
using DoseCalcProgressCallback = std::function<void(int, int, const std::string&)>;
// (currentBeam, totalBeams, message)

class IDoseEngine {
    virtual std::string getName() const = 0;

    virtual DoseInfluenceMatrix calculateDij(
        const Plan& plan, const Stf& stf,
        const PatientData& patientData, const Grid& doseGrid) = 0;

    virtual DoseMatrix calculateDose(
        const DoseInfluenceMatrix& dij,
        const std::vector<double>& weights,
        const Grid& grid) = 0;

    void setProgressCallback(DoseCalcProgressCallback cb);
    void setCancelFlag(std::atomic<bool>* flag);
    void setOptions(const DoseCalcOptions& opts);
};
```

### DoseCalcOptions

```cpp
struct DoseCalcOptions {
    double absoluteThreshold = 1e-6;    // Minimum absolute dose to store (Gy)
    double relativeThreshold = 0.01;    // Minimum relative dose fraction per bixel
    int numThreads = 0;                 // 0 = use all available threads
};
```

## Dose Engine

### PencilBeamEngine

Implements the Bortfeld (1993) SVD pencil beam algorithm with FFT-convolved lateral kernels.

**Inherits:** `IDoseEngine`

```cpp
class PencilBeamEngine : public IDoseEngine {
    std::string getName() const override;   // "PencilBeam"

    DoseInfluenceMatrix calculateDij(
        const Plan& plan, const Stf& stf,
        const PatientData& patientData, const Grid& doseGrid) override;

    DoseMatrix calculateDose(
        const DoseInfluenceMatrix& dij,
        const std::vector<double>& weights,
        const Grid& grid) override;
};
```

#### Algorithm Overview

The pencil beam algorithm decomposes dose into depth and lateral components:

$$D(x, y, z) = \sum_{i=1}^{3} w_i \cdot D_{\text{depth},i}(z) \cdot K_{\text{lateral},i}(x, y, z)$$

Where the three kernel components are:
1. **Primary** ($\beta_1$) — Direct photon interactions
2. **First scatter** ($\beta_2$) — Single-scattered photons
3. **Second scatter** ($\beta_3$) — Multiple-scattered photons

**Steps per beam:**

1. **Ray tracing** — Use `SiddonRayTracer` to find voxel intersections along each ray
2. **SSD computation** — Use `SSDCalculator` for source-to-surface distance per ray
3. **Radiological depth** — Use `RadDepthCalculator` to compute water-equivalent depths through electron density
4. **Depth dose** — Apply Bortfeld analytical depth-dose model with SVD decomposition
5. **Lateral convolution** — Convolve fluence with lateral spread kernels using `FFT2D`
6. **Assembly** — Populate `DoseInfluenceMatrix` with (voxel, bixel, dose) triplets

#### Machine Parameters Used

| Parameter | Source | Description |
|-----------|--------|-------------|
| `betas[3]` | `MachineData` | Lateral kernel widths (σ) |
| `m` | `MachineData` | Scatter-to-primary ratio |
| `primaryFluence` | `MachineData` | Off-axis fluence profile |
| `kernel[].kernelN` | `MachineData` | Depth-dependent kernel profiles |
| `penumbraFWHMatIso` | `MachineData` | Penumbra full width at half max |
| `SAD`, `SCD` | `MachineMeta` | Source distances |

### DoseEngineFactory

```cpp
class DoseEngineFactory {
    static std::unique_ptr<IDoseEngine> create(const std::string& engineName);
};
```

Currently supports: `"PencilBeam"`.

## Dose Influence Matrix

### DoseInfluenceMatrix

Sparse matrix mapping bixel weights to voxel doses: $\vec{d} = D_{ij} \cdot \vec{w}$

Uses a two-phase construction pattern:

```cpp
class DoseInfluenceMatrix {
    DoseInfluenceMatrix(size_t numVoxels, size_t numBixels);

    // Phase 1: COO construction
    void setValue(size_t voxel, size_t bixel, double value);
    void reserveNonZeros(size_t nnz);
    void appendBatch(const std::vector<size_t>& rows,
                     const std::vector<size_t>& cols,
                     const std::vector<double>& vals);

    // Phase transition
    void finalize();                   // COO → CSR conversion
    bool isFinalized() const;

    // Phase 2: CSR access
    double getValue(size_t voxel, size_t bixel) const;
    double operator()(size_t voxel, size_t bixel) const;

    // Key operations
    std::vector<double> computeDose(const std::vector<double>& weights) const;
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                    std::vector<double>& grad) const;

    // Dimensions
    size_t getNumVoxels() const;
    size_t getNumBixels() const;
    size_t getNumNonZeros() const;

    // Direct CSR data access (for serialization)
    const std::vector<double>&  getValues() const;
    const std::vector<size_t>&  getColIndices() const;
    const std::vector<size_t>&  getRowPtrs() const;
    void loadCSR(std::vector<size_t> rowPtrs,
                 std::vector<size_t> colIndices,
                 std::vector<double>  values);
};
```

#### COO → CSR Storage

**COO (Coordinate) format** — Used during construction:
- Triplets `(row, col, value)` appended in arbitrary order
- Thread-safe batch append via `appendBatch()`

**CSR (Compressed Sparse Row) format** — Used for computation:
- `rowPtrs[i]` = start of row `i` in `colIndices`/`values`
- Efficient row-wise access for `computeDose()` and transpose product

#### Key Operations

**Forward:** Compute dose from weights (used in optimization loop):
```cpp
std::vector<double> dose = dij.computeDose(weights);
// dose[v] = Σ_b dij(v,b) × weights[b]
```

**Transpose:** Compute gradient w.r.t. weights (used for L-BFGS-B):
```cpp
dij.accumulateTransposeProduct(voxelGradient, weightGradient);
// weightGradient[b] += Σ_v dij(v,b) × voxelGradient[v]
```

### DijSerializer

Binary serialization with deterministic caching.

```cpp
class DijSerializer {
    static bool save(const DoseInfluenceMatrix& dij, const std::string& filePath);
    static DoseInfluenceMatrix load(const std::string& filePath);
    static bool exists(const std::string& filePath);

    static std::string buildCacheKey(const std::string& patientName,
                                      int numBeams, double bixelWidth, double doseResX);
    static std::string getCacheDir();
};
```

**Binary format (v2):**
- Magic: `"ODIJ"` (4 bytes)
- Version: `2` (uint32)
- Dimensions: numVoxels, numBixels, numNonZeros (uint64 each)
- CSR arrays: rowPtrs, colIndices, values (raw binary)

**Cache key format:** `<PatientName>_<numBeams>beams_bw<bixelWidth>_res<resolution>mm.dij`

Example: `DOE^JOHN_90beams_bw7.0_res2.5mm.dij`

### DoseMatrix

3D dose result on a grid.

```cpp
class DoseMatrix {
    void setGrid(const Grid& grid);
    void allocate();
    double& at(size_t i, size_t j, size_t k);
    double getMax() const;
    double getMean() const;
    double* data();
    size_t size() const;
};
```

## Ray Tracing

### SiddonRayTracer

Implementation of the Siddon (1985) parametric ray-voxel intersection algorithm.

```cpp
struct RayTraceResult {
    std::vector<size_t> voxelIndices;
    std::vector<double> intersectionLengths;   // mm
    std::vector<double> densities;
    double totalDistance = 0.0;                 // mm
};

class SiddonRayTracer {
    static RayTraceResult trace(
        const Vec3& source, const Vec3& target,
        const Grid& grid, const double* densityData = nullptr);

    static std::vector<std::pair<size_t, double>> traceRadDepth(
        const Vec3& source, const Vec3& target,
        const Grid& grid, const double* densityData);
};
```

**Algorithm:**
1. Compute parametric intersections of the ray with all grid planes (X, Y, Z)
2. Merge intersection parameters into sorted sequence
3. For each segment between consecutive intersections, compute the voxel index and intersection length
4. Optionally accumulate radiological depth (length × density)

### SSDCalculator

Source-to-surface distance computation.

```cpp
class SSDCalculator {
    static double computeSSD(
        const Vec3& source, const Vec3& target,
        const Grid& grid, const double* densityData,
        double densityThreshold = 0.05);

    static std::vector<double> computeBeamSSDs(
        const Vec3& sourcePoint,
        const std::vector<Vec3>& rayTargets,
        const Grid& grid, const double* densityData,
        double densityThreshold = 0.05);
};
```

The SSD is the distance from the source to the first voxel where electron density exceeds the threshold (patient surface entry point).

### RadDepthCalculator

Radiological depth (water-equivalent path length) computation.

```cpp
class RadDepthCalculator {
    static std::unordered_map<size_t, double> computeRadDepths(
        const Vec3& sourcePoint, const Grid& grid,
        const double* densityData,
        const std::vector<size_t>& voxelIndices);

    static std::vector<std::pair<size_t, double>> computeRayRadDepths(
        const Vec3& sourcePoint, const Vec3& targetPoint,
        const Grid& grid, const double* densityData);
};
```

**Radiological depth** at a voxel is the integral of electron density along the ray from the surface to that voxel:

$$z_{\text{rad}} = \int_{\text{surface}}^{z} \rho_e(s) \, ds$$

## Plan Analysis

### PlanAnalysis

Per-structure dose statistics and DVH computation.

```cpp
struct StructureDoseStats {
    std::string name, type;
    size_t numVoxels;
    double minDose, maxDose, meanDose, stdDose;
    double d2, d5, d50, d95, d98;          // Dx% (Gy)
    double v20, v40, v50, v60;             // VxGy (%)
    double conformityIndex, homogeneityIndex;
    std::vector<double> sortedDoses;
};

struct DVHCurveData {
    std::string structureName;
    std::vector<float> doses;              // Gy
    std::vector<float> volumes;            // %
};

class PlanAnalysis {
    static std::vector<StructureDoseStats> computeStats(
        const DoseMatrix& dose, const PatientData& patient,
        const Grid& doseGrid, double prescribedDose = 60.0);

    static std::vector<DVHCurveData> computeDVHCurves(
        const std::vector<StructureDoseStats>& stats,
        double maxDose, int numBins = 200);

    static void print(const std::vector<StructureDoseStats>& stats,
                      std::ostream& os = std::cout);

    static double computeDx(const std::vector<double>& sortedDoses, double percentile);
    static double computeVx(const std::vector<double>& sortedDoses, double doseThreshold);
};
```

#### Quality Indices

| Metric | Definition | Description |
|--------|-----------|-------------|
| **Dx%** | Minimum dose to x% of the volume | D95 = dose received by 95% of the volume |
| **VxGy** | Percentage of volume receiving ≥ x Gy | V20 = % receiving ≥ 20 Gy |
| **CI** | $\frac{V_{\text{target receiving } \geq \text{Rx}}}{V_{\text{target}}}$ | Conformity Index (ideal = 1.0) |
| **HI** | $\frac{D_2 - D_{98}}{D_{50}}$ | Homogeneity Index (ideal = 0.0) |

## Mathematical Tools

### FFT2D

Cooley-Tukey 2D FFT for lateral kernel convolution.

```cpp
class FFT2D {
    using Complex = std::complex<double>;

    static void fft2(std::vector<Complex>& data, size_t rows, size_t cols, bool inverse = false);

    static std::vector<double> convolve2D(
        const std::vector<double>& a, size_t aRows, size_t aCols,
        const std::vector<double>& b, size_t bRows, size_t bCols);

    static std::vector<double> convolve2DSame(
        const std::vector<double>& a, size_t aRows, size_t aCols,
        const std::vector<double>& b, size_t bRows, size_t bCols);

    static size_t nextPow2(size_t n);
};
```

`convolve2DSame()` returns a result with the same dimensions as input `a` (centered crop of the full convolution).

### GridInterpolant2D

Bilinear interpolation on a regular 2D grid.

```cpp
class GridInterpolant2D {
    GridInterpolant2D(double xMin, double xMax, size_t nx,
                      double yMin, double yMax, size_t ny,
                      const std::vector<double>& data);

    double operator()(double x, double y) const;   // Returns 0 outside bounds

    std::vector<double> evaluate(const std::vector<double>& xs,
                                  const std::vector<double>& ys) const;
    bool isValid() const;
};
```

Used by the pencil beam engine for interpolating depth-dependent kernel values at arbitrary positions.

## Usage Example

```cpp
// Create engine
auto engine = DoseEngineFactory::create("PencilBeam");
engine->setOptions({.absoluteThreshold = 1e-6, .relativeThreshold = 0.01});
engine->setProgressCallback([](int beam, int total, const std::string& msg) {
    std::cout << "Beam " << beam << "/" << total << ": " << msg << "\n";
});

// Calculate Dij
auto doseGrid = Grid::createDoseGrid(patientData->getGrid(), {2.5, 2.5, 2.5});
auto dij = engine->calculateDij(*plan, *stf, *patientData, doseGrid);

// Check cache
std::string cacheKey = DijSerializer::buildCacheKey("DOE^JOHN", 9, 7.0, 2.5);
DijSerializer::save(dij, cacheKey);

// Compute dose from optimized weights
auto doseMatrix = engine->calculateDose(dij, weights, doseGrid);

// Analyze
auto stats = PlanAnalysis::computeStats(doseMatrix, *patientData, doseGrid, 60.0);
PlanAnalysis::print(stats);
auto dvhCurves = PlanAnalysis::computeDVHCurves(stats, doseMatrix.getMax());
```

## Related Documentation

- [Architecture](architecture.md) — Pipeline overview and module dependencies
- [Optimization Module](optimization.md) — Uses Dij for fluence optimization
- [Geometry Module](geometry.md) — Grid and Volume used for dose calculation
