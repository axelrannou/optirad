# Core Module (`optirad_core`)

The core module defines the central domain model for OptiRad: patient data, treatment plans, beams, machines, steering files, fluence containers, and workflow-facing sequencing types.

**Library:** `optirad_core`  
**Dependencies:** `optirad_utils`, `optirad_geometry`  
**Optional:** OpenMP (parallel STF generation)

## Files

| File | Description |
|------|-------------|
| `Patient.hpp/cpp` | Patient demographics |
| `PatientData.hpp/cpp` | Central data container |
| `Plan.hpp/cpp` | Treatment plan |
| `Beam.hpp/cpp` | Beam definition and ray generation |
| `Ray.hpp` | Single beamlet/ray |
| `Aperture.hpp` | Step-and-shoot aperture and sequencing result types |
| `FluenceMap.hpp` | 2D fluence map extracted from optimizer weights |
| `Machine.hpp/cpp` | Linear accelerator parameters |
| `Stf.hpp/cpp` | Steering file (beam collection) |
| `workflow/*.hpp/cpp` | Shared orchestration pipelines (`optirad_workflow`) |

## Classes

### Patient

Simple value type holding patient demographics.

```cpp
class Patient {
    void setName(const std::string& name);
    void setID(const std::string& id);
    const std::string& getName() const;
    const std::string& getID() const;
};
```

### PatientData

Central container for all patient-related data. Analogous to matRad's `ct` + `cst` structure.

```cpp
class PatientData {
    // Patient info
    void setPatient(std::unique_ptr<Patient> patient);
    Patient* getPatient();

    // CT volume (HU values, int16_t)
    void setCTVolume(std::unique_ptr<Volume<int16_t>> ct);
    Volume<int16_t>* getCTVolume();

    // Electron density volume (for dose calculation)
    void setEDVolume(std::unique_ptr<Volume<double>> ed);
    Volume<double>* getEDVolume();

    // Structures (ROIs)
    void setStructureSet(std::unique_ptr<StructureSet> structures);
    StructureSet* getStructureSet();

    // Grid from CT
    const Grid& getGrid() const;

    // Validation
    bool hasValidCT() const;
    bool hasStructures() const;
    bool isValid() const;

    // HU → electron density conversion (piecewise linear HLUT)
    void convertHUtoED();
};
```

**Ownership model:** `PatientData` owns `Patient`, CT `Volume`, ED `Volume`, and `StructureSet` via `std::unique_ptr`. It is shared across the pipeline via `std::shared_ptr<PatientData>`.

#### HU → Electron Density Conversion

The `convertHUtoED()` method applies a piecewise linear HLUT (Hounsfield Lookup Table):

| HU Range | Electron Density | Tissue |
|----------|-----------------|--------|
| HU ≤ -1000 | ED = 0.0 | Air |
| -1000 < HU ≤ 0 | ED = 1.0 + HU/1000 | Lung / soft tissue |
| 0 < HU ≤ 100 | ED = 1.0 + 0.001 × HU | Soft tissue |
| HU > 100 | ED = 1.1 + 0.0005 × (HU - 100) | Bone |

### Plan

Defines a treatment plan configuration.

```cpp
class Plan {
    void setName(const std::string& name);
    void setRadiationMode(const std::string& mode);   // "photons" (default)
    void setMachine(const Machine& machine);
    void setNumOfFractions(int n);                     // default: 30
    void setStfProperties(const StfProperties& stf);
    void setPatientData(std::shared_ptr<PatientData> data);

    void addBeam(const Beam& beam);
    const std::vector<Beam>& getBeams() const;
    size_t getNumBeams() const;

    // Compute isocenter from target structure center of mass
    std::array<double, 3> computeIsoCenter() const;

    void setDoseGridResolution(const std::array<double, 3>& res);  // default: {2.5, 2.5, 2.5} mm
    void printSummary() const;
};
```

### Beam

Represents a single treatment beam with angular geometry and a collection of rays/bixels.

```cpp
class Beam {
    // Angular parameters
    void setGantryAngle(double angle);          // degrees
    void setCouchAngle(double angle);           // degrees

    // Geometry
    void setIsocenter(const Vec3& iso);         // LPS coordinates (mm)
    void setBixelWidth(double width);           // mm (default: 7.0)
    void setSAD(double sad);                    // mm (default: 1000.0)
    void setSCD(double scd);                    // mm (default: 500.0)

    // Source points
    const Vec3& getSourcePointBev() const;      // Source in BEV coordinates
    const Vec3& getSourcePoint() const;         // Source in LPS coordinates
    void computeSourcePoints();                 // Compute from gantry/couch angles

    // Rays
    void addRay(const Ray& ray);
    const std::vector<Ray>& getRays() const;
    size_t getNumOfRays() const;
    size_t getTotalNumOfBixels() const;

    // Ray generation
    void generateRays(double bixelWidth, const std::array<double, 2>& fieldSize);
    void generateRaysFromTarget(const std::vector<Vec3>& targetWorldCoords,
                                double bixelWidth,
                                const Vec3& ctResolution);
    void computePhotonRayCorners();             // Compute beamlet corners at iso and SCD planes
    void setAllRayEnergies(double energy);
};
```

#### Beam Geometry

- **SAD** (Source-to-Axis Distance): Distance from radiation source to isocenter (typically 1000 mm)
- **SCD** (Source-to-Collimator Distance): Distance from source to collimator plane (typically 500 mm)
- **BEV** (Beam's Eye View): Coordinate system looking from source through isocenter
- **LPS** (Left-Posterior-Superior): DICOM patient coordinate system

The `computeSourcePoints()` method uses the gantry and couch rotation matrix (from `geometry/MathUtils.hpp`) to transform the source position from BEV to LPS coordinates.

### Ray

Single beamlet within a beam. Each ray defines a pencil beam direction.

```cpp
class Ray {
    // BEV coordinates
    void setRayPosBev(const Vec3& pos);
    void setTargetPointBev(const Vec3& target);

    // LPS (patient) coordinates
    void setRayPos(const Vec3& pos);
    void setTargetPoint(const Vec3& target);

    // Beamlet geometry (photon-specific)
    void setBeamletCornersAtIso(const std::array<Vec3, 4>& corners);
    void setRayCornersSCD(const std::array<Vec3, 4>& corners);

    // Energy
    void setEnergy(double energy);              // default: 6.0 MV
    double getEnergy() const;

    size_t getNumOfBixels() const;              // Always returns 1 for photons
};
```

### Aperture and Sequencing Types

The core module also defines the data structures shared between optimization, sequencing, and GUI BEV rendering.

```cpp
struct Aperture {
    std::vector<double> bankA;
    std::vector<double> bankB;
    double weight = 0.0;
    int segmentIndex = 0;
};

struct LeafSequenceResult {
    std::vector<Aperture> segments;
    size_t beamIndex = 0;
    double totalMU = 0.0;
    double fluenceFidelity = 0.0;
    std::vector<int> quantizedFluence;
    double calFac = 0.0;
    int numLevels = 0;
    std::vector<double> leafPairBoundariesZ;
    std::vector<double> leafPairFluence;
    int leafPairFluenceCols = 0;
    double originX = 0.0;
};

struct LeafSequencerOptions {
    int numLevels = 15;
    double minSegmentMU = 0.0;
    double leafPositionResolution = 0.5;
};
```

- `Aperture` stores one step-and-shoot segment as a pair of leaf-bank openings at isocenter plane plus MU weight.
- `LeafSequenceResult` stores the per-beam deliverable representation, fidelity metrics, and fluence data reused by the sequencing module and `BevView`.
- `LeafSequencerOptions` configures intensity quantization, minimum segment MU filtering, and leaf-position snapping.

### FluenceMap

`FluenceMap` is a lightweight 2D BEV grid extracted from the global optimizer weight vector for one beam.

```cpp
class FluenceMap {
public:
    static FluenceMap fromBeamWeights(const Beam& beam,
                                      const std::vector<double>& weights,
                                      size_t globalOffset);

    double getValue(int row, int col) const;
    std::vector<double> getProfile(int row) const;
    std::vector<int> mapToLeafPairs(const MachineGeometry& mlc) const;
    double getMaxFluence() const;

    int getNumRows() const;
    int getNumCols() const;
    double getBixelWidth() const;
    double getOriginX() const;
    double getOriginZ() const;
};
```

The map is organized in BEV coordinates:

- rows correspond to beam-space `Z`,
- columns correspond to beam-space `X`,
- and `mapToLeafPairs()` bridges the optimizer grid to physical MLC leaf-pair geometry.

### Machine

Defines a linear accelerator's physical parameters. Two machine types are supported:

#### MachineType Enum

```cpp
enum class MachineType { Generic, PhaseSpace };
```

#### Machine Components

| Component | Description |
|-----------|-------------|
| `MachineMeta` | Basic metadata: name, radiation mode, SAD, SCD, machine type |
| `MachineData` | Physical parameters: beta values, energy, primary fluence profile, dose kernels, penumbra FWHM |
| `MachineConstraints` | Mechanical limits: gantry rotation speed, leaf speed, monitor unit rate |
| `MachineGeometry` | Geometric parameters: jaw ranges, MLC configuration, field sizes, phase-space file paths |

```cpp
class Machine {
    const MachineMeta& getMeta() const;
    const MachineData& getData() const;
    const MachineConstraints& getConstraints() const;
    const MachineGeometry& getGeometry() const;

    // Convenience
    const std::string& getName() const;
    double getSAD() const;
    double getSCD() const;
    MachineType getMachineType() const;
    bool isPhaseSpace() const;
};
```

#### Supported Machines

| Machine | Type | Configuration |
|---------|------|--------------|
| Generic | `MachineType::Generic` | Single JSON file (`machine_photons_Generic.json`) |
| Varian_TrueBeam6MV | `MachineType::PhaseSpace` | JSON config + 6 IAEA phase-space file pairs |

#### Kernel Data (MachineData)

The pencil beam engine uses three kernel components from `MachineData`:

- `betas[0]` (β₁) — Primary kernel width parameter
- `betas[1]` (β₂) — First scatter kernel width parameter
- `betas[2]` (β₃) — Second scatter kernel width parameter
- `m` — Scatter-to-primary ratio
- `primaryFluence` — Radial primary fluence profile
- `kernel[]` — Depth-dependent kernel profiles at various SSDs

### Stf (Steering File)

Manages the collection of beams for a treatment plan.

```cpp
class Stf {
    void addBeam(const Beam& beam);
    const Beam* getBeam(size_t index) const;
    size_t getCount() const;
    bool isEmpty() const;
    const std::vector<Beam>& getBeams() const;

    // Aggregates
    size_t getTotalNumOfRays() const;
    size_t getTotalNumOfBixels() const;

    // Batch operations
    void computeAllSourcePoints();
    void generateAllRays(double bixelWidth, const std::array<double, 2>& fieldSize);
    void printSummary() const;
};
```

## Workflow Orchestration (`optirad_workflow`)

The `src/core/workflow/` subdirectory is built as a separate orchestration library named `optirad_workflow`. It extracts the common planning pipelines used by both the CLI and the GUI.

### WorkflowState

`WorkflowState` is the shared pipeline state container used to carry data across planning stages.

```cpp
struct WorkflowState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;
    std::shared_ptr<Stf> stf;
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;
    std::shared_ptr<DoseInfluenceMatrix> dij;
    std::shared_ptr<Grid> doseGrid;
    std::shared_ptr<Grid> computeGrid;
    std::vector<double> optimizedWeights;
    std::shared_ptr<DoseMatrix> doseResult;
    std::vector<LeafSequenceResult> leafSequences;
    std::vector<double> deliverableWeights;
    std::shared_ptr<DoseMatrix> deliverableDoseResult;
    std::vector<StructureDoseStats> deliverableStats;
    DoseManager doseManager;

    bool dicomLoaded() const;
    bool planCreated() const;
    bool stfGenerated() const;
    bool phaseSpaceLoaded() const;
    bool dijComputed() const;
    bool optimizationDone() const;
    bool leafSequenceDone() const;
    bool doseAvailable() const;

    void syncSelectedDose();
    void resetPlan();
    void resetStf();
    void resetPhaseSpace();
    void resetDij();
    void resetLeafSequence();
    void resetOptimization();
    void resetAllDoses();
};
```

### Pipeline Builders

The workflow layer groups the main shared execution paths:

| Type | Purpose |
|------|---------|
| `PlanBuilder` | Build a `Plan`, `StfProperties`, and generic-machine `Stf` from a `PlanConfig` |
| `DoseCalculationPipeline` | Create dose grid, check cache, compute Dij, and save cache output |
| `OptimizationPipeline` | Build objectives, run the optimizer, compute forward dose, and generate statistics |
| `PhaseSpaceBuilder` | Build one `PhaseSpaceBeamSource` per beam for phase-space machines |
| `LeafSequencingPipeline` | Convert optimized weights into fluence maps, sequence apertures, compute deliverable dose, and report statistics |

Representative configuration types:

```cpp
struct PlanConfig {
    std::string radiationMode = "photons";
    std::string machineName = "Generic";
    int numFractions = 30;
    double bixelWidth = 7.0;
    std::vector<double> gantryAngles;
    std::vector<double> couchAngles;
};

struct DoseCalcPipelineOptions {
    std::array<double, 3> resolution = {2.5, 2.5, 2.5};
    bool useCache = true;
    double absoluteThreshold = 1e-6;
    double relativeThreshold = 0.01;
    int numThreads = 0;
};

struct OptimizationConfig {
    int maxIterations = 400;
    double tolerance = 1e-5;
    double targetDose = 66.0;
    bool ntoEnabled = true;
    double ntoThresholdPct = 1.04;
    double ntoPenalty = 2000.0;
    double spatialSmoothingWeight = 0.0;
    double l2RegWeight = 0.0;
    double l1RegWeight = 0.0;
};
```

These pipeline types keep CLI and GUI behavior aligned while isolating the treatment-planning business logic from presentation code.

## Usage Example

```cpp
// Create plan
auto plan = std::make_shared<Plan>();
plan->setRadiationMode("photons");
plan->setMachine(machine);
plan->setNumOfFractions(30);
plan->setPatientData(patientData);
plan->setDoseGridResolution({2.5, 2.5, 2.5});

// Configure STF properties
StfProperties stfProps;
stfProps.setGantryAngles(0, 40, 360);   // 9 beams every 40°
stfProps.bixelWidth = 7.0;
stfProps.setUniformIsoCenter(plan->computeIsoCenter());
plan->setStfProperties(stfProps);
```

## Related Documentation

- [Geometry Module](geometry.md) — Grid, Volume, coordinate systems used by core
- [Dose Module](dose.md) — Dij, dose storage, and `DoseManager` used by workflow state
- [Sequencing Module](sequencing.md) — Uses `Aperture`, `LeafSequenceResult`, and `FluenceMap`
- [Steering Module](steering.md) — STF generation from plan configuration
- [I/O Module](io.md) — Machine loading from JSON files
- [Architecture](architecture.md) — Global pipeline and module dependency graph
