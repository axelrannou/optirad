# Steering Module (`optirad_steering`)

The steering module handles generation of steering file (STF) properties and beam/ray configurations for treatment planning.

**Library:** `optirad_steering` (INTERFACE — header-only)  
**Dependencies:** `optirad_core`, `optirad_geometry` (at include level)

## Files

| File | Description |
|------|-------------|
| `StfProperties.hpp` | Lightweight STF configuration |
| `IStfGenerator.hpp` | Abstract STF generator interface |
| `PhotonIMRTStfGenerator.hpp` | Photon IMRT STF generator with target-aware rays |

## Classes

### StfProperties

Lightweight data structure holding beam angle configuration and bixel parameters.

```cpp
struct StfProperties {
    std::vector<double> gantryAngles;
    std::vector<double> couchAngles;
    size_t numOfBeams = 0;
    double bixelWidth = 7.0;                   // mm
    std::vector<std::array<double, 3>> isoCenters;

    // Helper: generate uniform gantry angles
    void setGantryAngles(double start, double step, double stopExclusive);

    // Helper: set same isocenter for all beams
    void setUniformIsoCenter(const std::array<double, 3>& iso);
};
```

**Example:** 9 equispaced beams:
```cpp
StfProperties props;
props.setGantryAngles(0, 40, 360);  // 0°, 40°, 80°, ..., 320°
props.bixelWidth = 7.0;
props.setUniformIsoCenter(isocenter);
// numOfBeams = 9
```

### IStfGenerator

Abstract interface for STF generation strategies.

```cpp
class IStfGenerator {
    virtual ~IStfGenerator() = default;
    virtual std::unique_ptr<StfProperties> generate() const = 0;
};
```

### PhotonIMRTStfGenerator

Photon IMRT steering file generator with target-aware ray placement.

**Inherits:** `IStfGenerator`

```cpp
class PhotonIMRTStfGenerator : public IStfGenerator {
    PhotonIMRTStfGenerator(double start = 0.0, double step = 60.0,
                           double stop = 360.0, double bixelWidth = 7.0,
                           const std::array<double, 3>& iso = {0, 0, 0});

    // Configuration
    void setMachine(const Machine& machine);
    void setFieldSize(const std::array<double, 2>& fieldSize);
    void setRadiationMode(const std::string& mode);
    void setTargetVoxelCoords(const std::vector<Vec3>& coords);
    void setCTResolution(const Vec3& res);
    void setGrid(const Grid& grid);
    void setStructureSet(const StructureSet& ss);

    // IStfGenerator override — basic properties only
    std::unique_ptr<StfProperties> generate() const override;

    // Full generation — beams with rays and geometry
    Stf generateStf() const;
};
```

#### Ray Generation Algorithm

The `generateStf()` method creates beams with target-aware ray positions:

1. **Collect target voxels** — Find all TARGET structures and collect their world-coordinate voxel centers
2. **Margin expansion** — Optionally apply 3D morphological dilation (`VoxelDilation`) to expand the target volume by a margin
3. **For each beam angle:**
   a. Create `Beam` with gantry/couch angles, isocenter, SAD/SCD from machine
   b. Compute source point via `computeSourcePoints()`
   c. Generate rays using `generateRaysFromTarget()`:
      - Project target voxels into BEV (Beam's Eye View)
      - Create a regular bixel grid covering the projected target extent
      - Each bixel center becomes a ray position
   d. Compute beamlet corners at isocenter and SCD planes
   e. Set ray energies from machine data
4. **Return** `Stf` containing all beams

#### Target-Aware vs Field-Size Ray Generation

| Method | Description | Use Case |
|--------|-------------|----------|
| `generateRaysFromTarget()` | Rays placed to cover projected target footprint | Default for IMRT — efficient bixel coverage |
| `generateRays()` | Regular grid over fixed field size | Fallback / uniform field |

**OpenMP parallelism:** Beam generation is parallelized across beams when OpenMP is available.

## Usage Example

```cpp
// Create generator
PhotonIMRTStfGenerator gen(0.0, 40.0, 360.0, 7.0, isocenter);
gen.setMachine(machine);
gen.setFieldSize({100.0, 100.0});
gen.setRadiationMode("photons");

// Set target voxels for target-aware ray generation
std::vector<Vec3> targetCoords;
for (auto idx : targetStructure->getVoxelIndices()) {
    targetCoords.push_back(grid.voxelToPatient(indexToVec3(idx)));
}
gen.setTargetVoxelCoords(targetCoords);
gen.setCTResolution(grid.getSpacing());

// Generate full STF
Stf stf = gen.generateStf();
stf.printSummary();
// Output: numBeams, total rays, total bixels per beam
```

## Integration with Pipeline

The steering module sits between plan creation and dose calculation:

```
Plan (with StfProperties)
    │
    ▼
PhotonIMRTStfGenerator::generateStf()
    │
    ▼
Stf (beams with rays)
    │
    ▼
PencilBeamEngine::calculateDij(plan, stf, ...)
```

For phase-space machines, the steering module is bypassed — `PhaseSpaceBeamSource` directly provides beam geometry.

## Related Documentation

- [Core Module](core.md) — STF, Beam, Ray, StfProperties data structures
- [Geometry Module](geometry.md) — VoxelDilation for target margin expansion
- [Dose Module](dose.md) — Uses Stf as input for dose calculation
