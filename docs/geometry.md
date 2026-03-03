# Geometry Module (`optirad_geometry`)

The geometry module provides spatial data structures and operations: 3D grids, volumes, structures with contour rasterization, coordinate systems, and mathematical utilities.

**Library:** `optirad_geometry`  
**Dependencies:** `optirad_utils`  
**Optional:** OpenMP (voxel dilation)

## Files

| File | Description |
|------|-------------|
| `Grid.hpp/cpp` | 3D grid with coordinate transforms |
| `Volume.hpp/cpp` | Templated 3D volume on a Grid |
| `Structure.hpp/cpp` | ROI with contours and voxel rasterization |
| `StructureSet.hpp/cpp` | Collection of structures |
| `CoordinateSystem.hpp/cpp` | Patient↔world, voxel↔world transforms |
| `MathUtils.hpp` | Vec3, Mat3, rotation matrices, BEV↔LPS transforms |
| `VoxelDilation.hpp` | 3D morphological dilation (header-only) |

## Classes

### Grid

Defines a 3D rectilinear grid with DICOM image orientation and coordinate transforms.

```cpp
class Grid {
    // Setup
    void setDimensions(size_t nx, size_t ny, size_t nz);
    void setSpacing(double dx, double dy, double dz);      // mm
    void setOrigin(const Vec3& origin);                     // mm, LPS
    void setPatientPosition(const std::string& pos);        // e.g. "HFS"
    void setImageOrientation(const std::array<double, 6>& orient);
    void setSliceThickness(double thickness);

    // Queries
    std::array<size_t, 3> getDimensions() const;            // {nx, ny, nz}
    Vec3 getSpacing() const;                                // {dx, dy, dz}
    Vec3 getOrigin() const;
    size_t getNumVoxels() const;                            // nx × ny × nz

    // Coordinate transforms
    Vec3 voxelToPatient(const Vec3& ijk) const;             // Voxel → LPS
    Vec3 patientToVoxel(const Vec3& lps) const;             // LPS → Voxel

    // Direction vectors (from image orientation)
    Vec3 getRowDirection() const;
    Vec3 getColumnDirection() const;
    Vec3 getSliceDirection() const;

    // Coordinate arrays for each axis
    std::vector<double> getXCoordinates() const;
    std::vector<double> getYCoordinates() const;
    std::vector<double> getZCoordinates() const;

    // Static helpers
    static Grid createDoseGrid(const Grid& ctGrid, const Vec3& doseResolution);
    static std::vector<size_t> mapVoxelIndices(const Grid& fromGrid, const Grid& toGrid,
                                                const std::vector<size_t>& fromIndices);
};
```

#### Coordinate System

OptiRad uses the DICOM LPS (Left-Posterior-Superior) patient coordinate system:

**Voxel-to-patient transform:**

The `voxelToPatient()` method applies the DICOM affine transform:

$$P = O + i \cdot \Delta x \cdot \hat{r} + j \cdot \Delta y \cdot \hat{c} + k \cdot \Delta z \cdot \hat{s}$$

Where:
- $O$ = grid origin (image position of first slice)
- $\Delta x, \Delta y, \Delta z$ = voxel spacing
- $\hat{r}, \hat{c}, \hat{s}$ = row, column, slice direction vectors

#### Dose Grid Creation

`Grid::createDoseGrid()` creates a dose calculation grid from the CT grid with different resolution:

```cpp
Grid doseGrid = Grid::createDoseGrid(ctGrid, {2.5, 2.5, 2.5});  // 2.5mm dose grid
```

Voxel indices can be mapped between grids:
```cpp
auto doseIndices = Grid::mapVoxelIndices(ctGrid, doseGrid, ctIndices);
```

### Volume\<T\>

Template class for 3D volumetric data stored on a `Grid`.

```cpp
template<typename T>
class Volume {
    void setGrid(const Grid& grid);
    const Grid& getGrid() const;
    void allocate();                                // Allocate memory for grid dimensions

    T& at(size_t i, size_t j, size_t k);            // 3D access
    const T& at(size_t i, size_t j, size_t k) const;

    T* data();                                      // Raw data pointer
    const T* data() const;
    size_t size() const;                            // Total number of voxels

    std::vector<double> getXCoordinates() const;
    std::vector<double> getYCoordinates() const;
    std::vector<double> getZCoordinates() const;
};
```

**Type aliases:**
```cpp
using CTVolume   = Volume<int16_t>;    // Hounsfield Units
using DoseVolume = Volume<double>;      // Dose in Gy
```

**Memory layout:** Data is stored as a flat 1D array in row-major order (i varying fastest). Index: `k * ny * nx + j * nx + i`.

### Structure

Represents a single Region of Interest (ROI) with contour data and rasterized voxel indices.

```cpp
struct Contour {
    std::vector<std::array<double, 3>> points;     // LPS coordinates (mm)
    double zPosition = 0.0;                         // Slice position
};

class Structure {
    void setName(const std::string& name);
    void setType(const std::string& type);          // "TARGET", "OAR", "EXTERNAL", "UNKNOWN"
    void setROINumber(int num);
    void setColor(uint8_t r, uint8_t g, uint8_t b);

    // Contour data (from DICOM RT-STRUCT)
    void addContour(const Contour& contour);
    const std::vector<Contour>& getContours() const;
    size_t getContourCount() const;

    // Voxel indices (after rasterization)
    void setVoxelIndices(const std::vector<size_t>& indices);
    const std::vector<size_t>& getVoxelIndices() const;
    size_t getVoxelCount() const;

    // Optimization parameters
    void setPriority(int priority);
    void setAlphaX(double alpha);                   // Radiobiology parameter
    void setBetaX(double beta);                     // Radiobiology parameter

    // Visibility (GUI)
    void setVisible(bool visible);
    bool isVisible() const;

    // Rasterize contours to voxel indices
    void rasterizeContours(const Grid& ctGrid);
};
```

#### Structure Types

| Type | Description | Typical Use |
|------|-------------|-------------|
| `TARGET` | Target volume (PTV, CTV, GTV) | Optimization target with prescribed dose |
| `OAR` | Organ at risk | Dose constraint (overdose penalty) |
| `EXTERNAL` | External body contour | Patient outline |
| `UNKNOWN` | Unclassified structure | — |

#### Contour Rasterization

`rasterizeContours()` converts polygon contours to a set of voxel indices on the CT grid using a scanline algorithm:

1. For each contour slice, find the matching CT slice by z-position
2. Apply point-in-polygon test for each voxel center on that slice
3. Collect all enclosed voxel linear indices

### StructureSet

Collection of `Structure` objects with batch operations.

```cpp
class StructureSet {
    void addStructure(std::unique_ptr<Structure> structure);
    const Structure* getStructure(size_t index) const;
    Structure* getStructure(size_t index);
    const Structure* getStructureByName(const std::string& name) const;
    size_t getCount() const;
    void clear();

    // Rasterize all structure contours on the grid
    void rasterizeContours(const Grid& ctGrid);
};
```

**Ownership:** `StructureSet` owns `Structure` objects via `std::unique_ptr`. Move-only (copy deleted).

### CoordinateSystem

Patient↔world and voxel↔world coordinate transformations.

```cpp
class CoordinateSystem {
    Vec3 patientToWorld(const Vec3& patient) const;
    Vec3 worldToPatient(const Vec3& world) const;
    Vec3 voxelToWorld(int i, int j, int k) const;
    void worldToVoxel(const Vec3& world, int& i, int& j, int& k) const;

    void setOrigin(const Vec3& origin);
    void setSpacing(const Vec3& spacing);
};
```

### Mathematical Utilities (`geometry/MathUtils.hpp`)

#### Types

```cpp
using Vec3 = std::array<double, 3>;

struct Mat3 {
    double m[3][3];
    Mat3();                     // Identity matrix
    Vec3 operator*(const Vec3& v) const;    // Matrix-vector product
};
```

#### Vector Operations

```cpp
double dot(const Vec3& a, const Vec3& b);
Vec3   cross(const Vec3& a, const Vec3& b);
double norm(const Vec3& v);
Vec3   normalize(const Vec3& v);
Vec3   vecAdd(const Vec3& a, const Vec3& b);
Vec3   vecSub(const Vec3& a, const Vec3& b);
Vec3   vecScale(const Vec3& v, double s);
```

#### Matrix Operations

```cpp
Mat3 transpose(const Mat3& m);
Mat3 matMul(const Mat3& a, const Mat3& b);
Mat3 inverse(const Mat3& m);
```

#### Rotation Matrix

```cpp
Mat3 getRotationMatrix(double gantryAngleDeg, double couchAngleDeg);
```

Returns the combined gantry + couch rotation matrix used to transform between BEV (Beam's Eye View) and LPS (patient) coordinates.

#### BEV ↔ LPS Transforms

```cpp
Vec3 vecMulMatTranspose(const Vec3& v, const Mat3& m);
```

Used for transforming ray positions and directions between beam and patient coordinate systems.

### Voxel Dilation (`VoxelDilation.hpp`)

3D morphological dilation with 26-connectivity for margin expansion.

```cpp
std::vector<size_t> dilateVoxels(
    const std::vector<size_t>& targetVoxelIndices,     // 1-based indices
    const std::set<size_t>& allVoxelIndices,           // 1-based, patient surface
    const std::array<size_t, 3>& dimensions,           // {ny, nx, nz}
    const std::array<double, 3>& spacing,              // {dy, dx, dz} mm
    const std::array<double, 3>& margin);              // {mx, my, mz} mm
```

**Algorithm:**
1. Compute dilation steps per axis: `steps_i = ceil(margin_i / spacing_i)`
2. For each iteration, expand the voxel set by adding all 26-connected neighbors
3. Only include neighbors that are within the patient (in `allVoxelIndices`)
4. OpenMP-parallelized neighbor search

**Use case:** Expanding target volume for ray generation in `PhotonIMRTStfGenerator`.

## Related Documentation

- [Core Module](core.md) — Uses Grid/ Volume for PatientData, beams use MathUtils for BEV↔LPS
- [Dose Module](dose.md) — SiddonRayTracer operates on the Grid/Volume
- [Steering Module](steering.md) — VoxelDilation used for target margin expansion
