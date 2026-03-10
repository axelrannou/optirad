# I/O Module (`optirad_io`)

The I/O module handles data import and export: DICOM-RT file processing, RT-STRUCT parsing, and machine configuration loading.

**Library:** `optirad_io`  
**Dependencies:** `optirad_core`, `optirad_geometry`, `nlohmann_json`  
**Optional:** DCMTK (DICOM support), TBB (parallel RT-STRUCT parsing)

## Files

| File | Description |
|------|-------------|
| `IDataImporter.hpp` | Abstract importer interface |
| `IDataExporter.hpp` | Abstract exporter interface |
| `DicomImporter.hpp/cpp` | DICOM import (DCMTK-based) |
| `DicomExporter.hpp/cpp` | DICOM RT export |
| `RTStructParser.hpp/cpp` | RT-STRUCT file parser |
| `MachineLoader.hpp/cpp` | Machine JSON loader |

## Interfaces

### IDataImporter

Abstract interface for data import.

```cpp
class IDataImporter {
public:
    virtual ~IDataImporter() = default;
    virtual bool canImport(const std::string& path) const = 0;
    virtual std::unique_ptr<Patient> importPatient(const std::string& path) = 0;
    virtual std::unique_ptr<StructureSet> importStructures(const std::string& path) = 0;
};
```

### IDataExporter

Abstract interface for data export.

```cpp
class IDataExporter {
public:
    virtual ~IDataExporter() = default;
    virtual bool exportDose(const DoseVolume& dose, const std::string& path) = 0;
    virtual bool exportPlan(const Plan& plan, const std::string& path) = 0;
};
```

## Classes

### DicomImporter

DICOM-RT importer using the DCMTK library. Handles all four standard DICOM-RT object types.

**Inherits:** `IDataImporter`

```cpp
class DicomImporter : public IDataImporter {
    // IDataImporter overrides
    bool canImport(const std::string& path) const override;
    std::unique_ptr<Patient> importPatient(const std::string& path) override;
    std::unique_ptr<StructureSet> importStructures(const std::string& path) override;

    // Full import — loads everything from a DICOM directory
    std::unique_ptr<PatientData> importAll(const std::string& dirPath);

    // Individual loaders
    bool loadDirectory(const std::string& dirPath);
    bool loadCTSeries(const std::string& dirPath);
    bool loadRTStruct(const std::string& filePath);
    bool loadRTPlan(const std::string& filePath);
    bool loadRTDose(const std::string& filePath);

    std::unique_ptr<Volume<int16_t>> importCTVolume();
    std::unique_ptr<StructureSet> importStructuresWithContours();
};
```

#### Supported DICOM Types

| DICOM Type | SOP Class UID | Description |
|------------|---------------|-------------|
| CT Image Storage | `1.2.840.10008.5.1.4.1.1.2` | CT slices (HU values) |
| RT Structure Set | `1.2.840.10008.5.1.4.1.1.481.3` | ROI contours |
| RT Plan | `1.2.840.10008.5.1.4.1.1.481.5` | Treatment plan parameters |
| RT Dose | `1.2.840.10008.5.1.4.1.1.481.2` | Dose distribution |

#### Import Pipeline

```
DICOM Directory
    │
    ├── loadDirectory()          → Scan and classify DICOM files
    │
    ├── loadCTSeries()           → Sort slices by position, build 3D volume
    │   └── importCTVolume()     → Create Volume<int16_t> with Grid
    │
    ├── loadRTStruct()           → Parse RT-STRUCT file
    │   └── importStructuresWithContours() → Create StructureSet with contours
    │
    ├── loadRTPlan()             → Parse beam parameters (optional)
    │
    └── loadRTDose()             → Parse dose grid (optional)
```

The `importAll()` convenience method executes the full pipeline and returns a complete `PatientData` object with HU→ED conversion already applied.

#### Conditional Compilation

When DCMTK is not available (`OPTIRAD_HAS_DCMTK` not defined), the importer provides stub implementations that log an error and return nullptr.

### DicomExporter

Exports treatment plan and dose data to DICOM-RT format.

```cpp
class DicomExporter {
    bool exportRTPlan(const Plan& plan, const std::string& outputPath);
    bool exportRTDose(const DoseMatrix& dose, const std::string& outputPath);
};
```

### RTStructParser

Specialized parser for DICOM RT-STRUCT files. Extracts ROI definitions, contour data, and structure properties.

```cpp
class RTStructParser {
    std::unique_ptr<StructureSet> parse(const std::string& filePath);
};
```

**Internal operations:**
1. Extract ROI names and numbers from the RT-STRUCT
2. Parse contour sequences for each ROI (contour points in LPS coordinates)
3. Extract structure colors (RGB)
4. Determine structure type (TARGET/OAR/EXTERNAL) from DICOM RT ROI type or name heuristics

**TBB parallelism:** When TBB is available (`OPTIRAD_HAS_TBB`), structure parsing is parallelized across ROIs using `tbb::parallel_for`.

### MachineLoader

Loads machine configurations from JSON files. All methods are static.

```cpp
class MachineLoader {
    // Load by radiation mode + machine name (auto-detects path)
    static Machine load(const std::string& radiationMode,
                        const std::string& machineName,
                        const std::string& dataDir);

    // Load from explicit file path
    static Machine loadFromFile(const std::string& filePath);

    // Load using compiled-in default data dir (OPTIRAD_DATA_DIR)
    static Machine load(const std::string& radiationMode,
                        const std::string& machineName);
};
```

#### Machine Type Auto-Detection

The loader auto-detects the machine type based on the file system layout:

| Layout | Machine Type | Description |
|--------|-------------|-------------|
| `data/machines/machine_<mode>_<name>.json` | `Generic` | Single JSON file with pencil beam parameters |
| `data/machines/<name>/` (directory) | `PhaseSpace` | JSON config + IAEA phase-space files |

#### Generic Machine JSON Format

```json
{
  "meta": {
    "radiationMode": "photons",
    "SAD": 1000.0,
    "SCD": 500.0,
    "machineType": "Generic"
  },
  "data": {
    "energy": 6.0,
    "betas": [0.3252, 0.0160, 0.0051],
    "m": 0.0051,
    "penumbraFWHMatIso": 5.0,
    "primaryFluence": [[0, 1.0], [50, 0.98], ...],
    "kernelPos": [0, 5, 10, ...],
    "kernel": [
      { "SSD": 900, "kernel1": [...], "kernel2": [...], "kernel3": [...] }
    ]
  },
  "constraints": { ... },
  "geometry": { ... }
}
```

#### Phase-Space Machine Layout

```
data/machines/Varian_TrueBeam6MV/
├── machine.json                 # Machine configuration
├── beam1.IAEAheader             # IAEA header file
├── beam1.IAEAphsp               # IAEA binary phase-space data
├── beam2.IAEAheader
├── beam2.IAEAphsp
└── ...
```

## Usage Example

```cpp
// Import DICOM
DicomImporter importer;
if (importer.canImport("/path/to/dicom")) {
    auto patientData = importer.importAll("/path/to/dicom");
}

// Load machine
auto machine = MachineLoader::load("photons", "Generic");
auto psMachine = MachineLoader::load("photons", "Varian_TrueBeam6MV");
```

## Related Documentation

- [Core Module](core.md) — PatientData, Machine structures populated by IO
- [Phase-Space Module](phsp.md) — IAEA file format details for phase-space machines
