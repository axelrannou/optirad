# Core Module (`optirad_core`)

The core module defines the central domain model for OptiRad: patient data, treatment plans, beams, machines, and steering files.

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
| `Machine.hpp/cpp` | Linear accelerator parameters |
| `Stf.hpp/cpp` | Steering file (beam collection) |

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
- [Steering Module](steering.md) — STF generation from plan configuration
- [I/O Module](io.md) — Machine loading from JSON files
