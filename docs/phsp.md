# Phase-Space Module (`optirad_phsp`)

The phase-space module provides IAEA phase-space file support for Monte Carlo beam sources: parsing header files, reading binary particle data, filtering, sampling, and beam source construction.

**Library:** `optirad_phsp`  
**Dependencies:** `optirad_core`, `optirad_geometry`  
**Optional:** OpenMP (parallel multi-file reading)

## Files

| File | Description |
|------|-------------|
| `IAEAHeaderParser.hpp/cpp` | IAEA header file parser |
| `IAEAPhspReader.hpp/cpp` | IAEA binary phase-space reader |
| `PhaseSpaceData.hpp/cpp` | Particle container with filtering/sampling |
| `PhaseSpaceBeamSource.hpp/cpp` | Full beam source pipeline |

## IAEA Phase-Space Format

OptiRad implements the IAEA(NDS)-0484 phase-space format specification. Each beam source consists of a pair of files:

| File | Extension | Content |
|------|-----------|---------|
| Header | `.IAEAheader` | ASCII metadata: particle counts, energy stats, record structure |
| Data | `.IAEAphsp` | Binary particle records: type, position, direction, energy, weight |

### Header Information

```cpp
struct IAEAHeaderInfo {
    int iaeaIndex;
    std::string title;
    int fileType;

    // Record structure flags
    bool storesEnergy, storesX, storesY, storesZ;
    bool storesU, storesV, storesWeight;
    int numExtraFloats, numExtraLongs;
    double constantWeight;
    int recordLength;
    int byteOrder;                          // 1234 = little-endian

    // Particle counts
    int64_t origHistories, totalParticles;
    int64_t numPhotons, numElectrons, numPositrons;

    // Energy statistics
    double photonEnergyCutoff, particleEnergyCutoff;
    double photonMeanEnergy, electronMeanEnergy, positronMeanEnergy;

    // Geometry extents (cm in file, converted to mm internally)
    std::array<double, 2> xRange, yRange, zRange;

    // Metadata
    std::string mcCodeVersion;
    double nominalSSD;
    std::string coordinateSystemDesc;

    int computeRecordLength() const;
    bool needsByteSwap() const;
};
```

### IAEAHeaderParser

```cpp
class IAEAHeaderParser {
    static IAEAHeaderInfo parse(const std::string& headerFilePath);
};
```

Parses the ASCII `.IAEAheader` file and extracts all metadata fields including the record structure (which fields are stored per particle, their byte offsets, and total record length).

## Particle Data

### PhaseSpaceData

Container for particle data with filtering and sampling operations.

```cpp
enum class ParticleType : int8_t {
    Photon   =  1,
    Electron = -1,
    Positron =  2
};

struct Particle {
    ParticleType type = ParticleType::Photon;
    std::array<double, 3> position  = {0, 0, 0};    // mm
    std::array<double, 3> direction = {0, 0, 1};     // direction cosines (unit vector)
    double energy = 0.0;                              // MeV
    double weight = 1.0;                              // statistical weight
};

struct PhaseSpaceMetrics {
    int64_t totalCount, photonCount, electronCount, positronCount;
    double minEnergy, maxEnergy, meanEnergy;
    double angularSpreadU, angularSpreadV;
    std::array<double, 2> xRange, yRange, zRange;    // mm
};

class PhaseSpaceData {
    std::vector<Particle>& particles();
    const std::vector<Particle>& particles() const;

    void addParticle(const Particle& p);
    void reserve(size_t n);
    size_t size() const;
    bool empty() const;
    void clear();

    // Random sampling (with seed for reproducibility)
    PhaseSpaceData sample(size_t n, unsigned int seed = 42) const;

    // In-place filtering
    void filterByJaws(double jawX1, double jawX2,
                      double jawY1, double jawY2,
                      double SAD, double scoringPlaneZ);
    void filterByType(ParticleType type);

    // Statistics
    PhaseSpaceMetrics computeMetrics() const;
    int64_t countByType(ParticleType type) const;
};
```

#### Jaw Filtering

`filterByJaws()` removes particles outside the field aperture defined by jaw positions. Particles are projected from the scoring plane to the isocenter plane using the SAD, then checked against jaw boundaries.

### IAEAPhspReader

Binary phase-space file reader. All methods are static.

```cpp
class IAEAPhspReader {
    // Read all particles from a single file
    static PhaseSpaceData readAll(const std::string& phspFilePath,
                                   const IAEAHeaderInfo& header);

    // Read a subset (offset + count)
    static PhaseSpaceData readSubset(const std::string& phspFilePath,
                                      const IAEAHeaderInfo& header,
                                      int64_t offset, int64_t count);

    // Read from multiple files (OpenMP parallel)
    static PhaseSpaceData readMultiple(const std::vector<std::string>& basePaths,
                                        const IAEAHeaderInfo& header,
                                        int64_t maxTotal = 0);

    // Random sampling without reading entire file
    static PhaseSpaceData readSampled(const std::string& phspFilePath,
                                       const IAEAHeaderInfo& header,
                                       int64_t sampleSize);
};
```

**Binary record format:** Each particle is stored as a fixed-length binary record. The record structure is defined by the header flags:

```
[particleType: int8] [energy: float32] [x: float32] [y: float32] [z: float32]
[u: float32] [v: float32] [weight: float32] [extra floats...] [extra longs...]
```

Fields may be omitted based on header flags (e.g., if `storesWeight = false`, the weight field uses `constantWeight`).

## Beam Source

### PhaseSpaceBeamSource

Full pipeline for constructing a treatment beam from phase-space data.

```cpp
class PhaseSpaceBeamSource {
    void configure(const Machine& machine,
                   double gantryAngle, double collimatorAngle,
                   double couchAngle,
                   const std::array<double, 3>& isocenter);

    void build(int64_t maxParticles = 0,          // 0 = all
               int64_t vizSampleSize = 100000);    // for visualization

    bool isBuilt() const;

    const PhaseSpaceData& getData() const;              // Full data
    const PhaseSpaceData& getVisualizationSample() const;  // Sampled for rendering
    const PhaseSpaceMetrics& getMetrics() const;
    const std::array<double, 3>& getIsocenter() const;
    const std::array<double, 3>& getSourcePosition() const;
    double getGantryAngle() const;
    const IAEAHeaderInfo& getHeaderInfo() const;

    std::vector<std::pair<double, int64_t>> computeEnergyHistogram(int numBins = 50) const;
};
```

#### Build Pipeline

```
configure()                          build()
    │                                   │
    ├── Store machine params            ├── Parse IAEA header(s)
    ├── Store beam angles               ├── Read phase-space data (all files)
    ├── Store isocenter                 ├── Filter by jaw aperture
    └── Compute source position         ├── Transform to patient LPS coords
        (gantry rotation)               ├── Compute metrics
                                        ├── Create visualization sample
                                        └── Compute energy histogram
```

#### Coordinate Transform

After reading, particle positions and directions are transformed from the machine coordinate system (beam-centered) to the patient LPS coordinate system using the gantry/couch rotation matrix:

$$\vec{p}_\text{LPS} = R(\theta_g, \theta_c) \cdot \vec{p}_\text{machine} + \vec{p}_\text{isocenter}$$

## Usage Example

```cpp
// Load phase-space machine
auto machine = MachineLoader::load("photons", "Varian_TrueBeam6MV");

// Create beam source
auto source = std::make_shared<PhaseSpaceBeamSource>();
source->configure(machine, 0.0, 0.0, 0.0, isocenter);
source->build(1000000, 100000);  // 1M particles, 100K for visualization

// Access data
auto metrics = source->getMetrics();
std::cout << "Photons: " << metrics.photonCount << "\n";
std::cout << "Mean energy: " << metrics.meanEnergy << " MeV\n";

auto histogram = source->computeEnergyHistogram(50);
```

## Supported Machine

| Machine | Phase-Space Files | Description |
|---------|------------------|-------------|
| Varian TrueBeam 6MV | 6 pairs (`.IAEAheader` + `.IAEAphsp`) | 6 MV photon beam, multiple field configurations |

## Related Documentation

- [Core Module](core.md) — Machine structure with phase-space geometry fields
- [I/O Module](io.md) — MachineLoader auto-detects PhaseSpace machines
- [GUI Module](gui.md) — PhaseSpacePanel and PhaseSpaceRenderer for visualization
