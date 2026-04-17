# Segmentation Module (`optirad_segmentation`)

The segmentation module currently provides automatic generation of a BODY/EXTERNAL structure from CT data when no external contour is available in the imported RT-STRUCT.

**Library:** `optirad_segmentation`
**Dependencies:** `optirad_geometry`

## Files

| File | Description |
|------|-------------|
| `BodyContourGenerator.hpp/cpp` | Automatic BODY contour generation from CT volume |

## Classes

### BodyContourGenerator

Generates a `Structure` representing the external patient contour directly from the CT voxel data.

```cpp
class BodyContourGenerator {
public:
    static constexpr int16_t DEFAULT_HU_THRESHOLD = -300;

    static std::unique_ptr<Structure> generate(
        const Volume<int16_t>& ctVolume,
        int16_t huThreshold = DEFAULT_HU_THRESHOLD);
};
```

The generated structure:

- is named `BODY`,
- has type `EXTERNAL`,
- includes precomputed voxel indices for dose and optimization workflows,
- and includes per-slice contour polygons for GUI visualization.

## Algorithm

For each axial slice, the generator performs the following steps:

1. Threshold the CT slice at a configurable HU value.
2. Flood-fill from the slice borders to identify exterior air.
3. Fill internal cavities so the body mask becomes a solid envelope.
4. Keep only the largest connected component to suppress couch, table, and small artifacts.
5. Convert the binary mask into 1-based voxel indices using the same indexing convention as `Structure` rasterization.
6. Extract a contour polygon suitable for visualization in slice views.

## API Details

### `generate()`

```cpp
static std::unique_ptr<Structure> generate(
    const Volume<int16_t>& ctVolume,
    int16_t huThreshold = DEFAULT_HU_THRESHOLD);
```

- `ctVolume`: input CT data in Hounsfield units
- `huThreshold`: voxels above this threshold are classified as body; default is `-300`
- return value: a new `Structure` ready to be inserted into a `StructureSet`

## Internal Processing Stages

The implementation is organized around a small set of slice-level helpers:

- `createSliceMask()` builds the initial binary mask for one slice.
- `floodFillExterior()` marks connected exterior air regions.
- `keepLargestComponent()` removes disconnected non-patient components.
- `collectSliceVoxelIndices()` converts the final mask into 1-based voxel indices.
- `extractSliceContour()` builds a non-self-intersecting contour polygon from the mask.

## Notes

- The implementation operates slice-by-slice in the axial direction.
- The default threshold is tuned for a robust BODY envelope rather than fine anatomical segmentation.
- Internal low-density regions such as lungs or air cavities are intentionally filled so the resulting structure represents the external patient boundary.
- This module currently focuses on external contour generation only; it is not a general-purpose anatomical segmentation framework.

## Usage Example

```cpp
if (!patientData->getStructureSet()->hasStructure("BODY")) {
    auto body = BodyContourGenerator::generate(*patientData->getCTVolume());
    patientData->getStructureSet()->addStructure(std::move(body));
}
```
