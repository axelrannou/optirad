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

    // Gantry angles: start/step/stop or explicit list
    void setGantryAngles(double start, double step, double stopExclusive);
    void setGantryAngles(const std::vector<double>& angles);

    // Couch angles: start/step/stop, explicit list, or uniform
    void setCouchAngles(double start, double step, double stopExclusive);
    void setCouchAngles(const std::vector<double>& angles);
    void setUniformCouchAngle(double angle);

    // Ensure gantry/couch lists match in size (pads couch if needed)
    void ensureConsistentAngles();
    bool isValid() const;

    void setUniformIsoCenter(const std::array<double, 3>& iso);
};
```

**Two modes for combining angles:**

1. **Explicit lists** (paired 1:1, like matRad): `beam[i] = (gantryAngles[i], couchAngles[i])`.
   Both vectors must have the same length.
2. **Start/step/stop** (Cartesian product, multi-arc): `setCouchAngles(start, step, stop)` with step > 0 
   performs a Cartesian product against the existing gantry angles.
   Total beams = numGantry × numCouch.

**Example:** 9 equispaced beams with uniform couch angle:
```cpp
StfProperties props;
props.setGantryAngles(0, 40, 360);  // 0°, 40°, 80°, ..., 320°
props.setUniformCouchAngle(0.0);     // couch = 0° for all beams
props.bixelWidth = 7.0;
props.setUniformIsoCenter(isocenter);
// numOfBeams = 9
```

**Example:** 4 beams with explicit paired lists:
```cpp
StfProperties props;
props.setGantryAngles({0.0, 90.0, 180.0, 270.0});
props.setCouchAngles({0.0, 10.0, 20.0, 30.0});
props.bixelWidth = 7.0;
props.setUniformIsoCenter(isocenter);
// beam 0: gantry=0°, couch=0°
// beam 1: gantry=90°, couch=10°
// ...
```

**Example:** Multi-arc with Cartesian product (start/step/stop):
```cpp
props.setGantryAngles(0, 4, 360);       // 90 gantry entries
props.setCouchAngles(-5.0, 5.0, 10.0);  // 3 couch entries: -5, 0, 5
// Cartesian product: 90 × 3 = 270 beams
// Arc 1: gantry=[0,4,...,356] couch=-5
// Arc 2: gantry=[0,4,...,356] couch=0
// Arc 3: gantry=[0,4,...,356] couch=5
props.setUniformIsoCenter(isocenter);
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

    // Couch angle configuration (paired 1:1 with gantry)
    void setCouchAngles(double start, double step, double stop);
    void setCouchAngles(const std::vector<double>& angles);

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

// Optional: set couch angles (paired with gantry)
gen.setCouchAngles({0, 0, 10, 10, 20, 20, 30, 30, 0});  // explicit list
// Or: gen.setCouchAngles(0.0, 5.0, 45.0);                // start/step/stop

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
