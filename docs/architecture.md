# OptiRad — Architecture Overview

This document describes the global architecture of OptiRad, a modular C++17 radiotherapy treatment planning system (TPS).

## System Overview

OptiRad provides a complete treatment planning pipeline: DICOM import, structure management, dose calculation, fluence optimization, leaf sequencing, and plan analysis. The source tree is organized into 11 primary modules plus a shared workflow orchestration library, consumed by CLI and GUI front-ends.

## Module Dependency Graph

```
optirad_utils           ← leaf library, no dependencies
    ↑
optirad_geometry        ← Grid, Volume, Structure, coordinate math
    ↑
├── optirad_segmentation ← BODY contour generation from CT
│
optirad_core            ← Patient, Plan, Beam, Machine, Stf, Aperture, FluenceMap
    ↑
├── optirad_io          ← DICOM import/export, machine loading
├── optirad_phsp        ← IAEA phase-space beam sources
├── optirad_steering    ← STF generation (header-only INTERFACE library)
└── optirad_dose        ← Dose engines, Dij, ray tracing, plan analysis, dose manager
    ↑
    ├── optirad_optimization ← L-BFGS-B optimizer, objectives, constraints
    └── optirad_sequencing   ← Leaf sequencing and deliverable dose

optirad_workflow        ← Shared plan / dose / optimize / sequence pipelines
    ↑
optirad_gui             ← GUI application
```

Each module is built as a separate CMake library (`add_library`), except `optirad_steering` which is an INTERFACE (header-only) library. In addition, `optirad_workflow` provides shared orchestration logic reused by the front-ends.

## Treatment Planning Pipeline

The end-to-end treatment planning workflow follows a linear pipeline:

```
 ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌──────────┐    ┌────────────┐    ┌─────────┐
 │ Load DICOM  │───▶│ Create Plan │───▶│ Generate STF│───▶│ Calculate   │───▶│ Optimize │───▶│ Sequence   │───▶│ Analyze │
 │ (io)        │    │ (core)      │    │ (steering)  │    │ Dij (dose)  │    │ (optim.) │    │ Leaves     │    │ (dose)  │
 └─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘    └──────────┘    └────────────┘    └─────────┘
     │                  │                  │                   │                 │                  │                │
 PatientData         Plan +             Stf / PSF           DoseInfluence     Optimized        Apertures +      PlanAnalysis
 (CT, ED,            Machine            beams               Matrix (Dij)      weights          deliverable      (DVH, Dx%,
 structures,                                                                                     dose             CI, HI)
 optional BODY)
```

**Step details:**

1. **Load DICOM** — `DicomImporter::importAll()` reads CT series, RT-STRUCT, RT-Plan, RT-Dose. HU values are converted to electron density, and `BodyContourGenerator` can create a BODY structure when no external contour is present.
2. **Create Plan** — `PlanBuilder` configures `Plan` with machine, gantry angles, bixel width, fractions, and isocenter.
3. **Generate STF / Build PSF** — `PhotonIMRTStfGenerator::generateStf()` creates beams with target-aware ray positions for generic machines, while `PhaseSpaceBuilder` creates one `PhaseSpaceBeamSource` per beam for phase-space machines.
4. **Calculate Dij** — `DoseCalculationPipeline` creates the dose grid, checks cache state, and calls `PencilBeamEngine::calculateDij()`. Results are cached via `DijSerializer`.
5. **Optimize** — `OptimizationPipeline` builds objectives, runs `LBFGSOptimizer`, computes forward dose, and stores statistics.
6. **Sequence Leaves** — `LeafSequencingPipeline` converts optimized weights into `FluenceMap`s, sequences apertures, and recomputes deliverable dose.
7. **Analyze** — `PlanAnalysis::computeStats()` and `computeDVHCurves()` generate per-structure statistics and DVH data for both optimized and deliverable dose maps.

## Design Patterns

### Strategy Pattern (Abstract Interfaces)

All major computation components are defined via abstract interfaces, enabling extension without modification:

| Interface | Purpose | Implementations |
|-----------|---------|-----------------|
| `IDoseEngine` | Dose calculation | `PencilBeamEngine` |
| `IOptimizer` | Fluence optimization | `LBFGSOptimizer` |
| `IDataImporter` | Data import | `DicomImporter` |
| `IDataExporter` | Data export | `DicomExporter` |
| `IStfGenerator` | STF generation | `PhotonIMRTStfGenerator` |
| `IPanel` | GUI panel | 10 panel implementations |
| `IView` | GUI view | `SliceView`, `View3D`, `DVHView`, `BevView` |

### Factory Pattern

Engine and optimizer creation is abstracted via factories:

```cpp
auto engine    = DoseEngineFactory::create("PencilBeam");
auto optimizer = OptimizerFactory::create("LBFGS");
```

### Template Method

`ObjectiveFunction` defines the abstract `compute()` and `gradient()` methods, with shared logic for weight, structure binding, and voxel index mapping. Four concrete implementations specialize the computation.

### COO → CSR Builder (DoseInfluenceMatrix)

The dose influence matrix uses a two-phase construction:
1. **Build phase (COO):** Threads call `appendBatch()` with (row, col, value) triplets — lock-free accumulation.
2. **Finalize phase (CSR):** `finalize()` converts to compressed sparse row format for efficient matrix-vector products.

### Pipeline Pattern

Both CLI and GUI orchestrate the planning workflow as a sequential pipeline, with each step producing output consumed by the next. The shared `optirad_workflow` library centralizes that business logic so the front-ends do not need to duplicate plan-building, dose calculation, optimization, phase-space loading, or sequencing behavior.

### Shared State

- **Workflow layer:** `WorkflowState` holds the shared planning data (`PatientData`, `Plan`, `Stf`, `Dij`, optimized weights, deliverable dose, `DoseManager`).
- **GUI:** `GuiAppState` extends that model with async task state, dose-selection caches, and view-refresh bookkeeping.

### PIMPL (Pointer to Implementation)

`View3D` uses the PIMPL idiom to hide OpenGL implementation details from the header.

## Build System

### CMake Configuration

- **Standard:** C++17 (`CMAKE_CXX_STANDARD 17`)
- **Minimum CMake:** 3.15
- **Compile commands:** Always exported (`CMAKE_EXPORT_COMPILE_COMMANDS ON`)

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `OPTIRAD_BUILD_TESTS` | `ON` | Build GoogleTest unit tests |
| `OPTIRAD_BUILD_GUI` | `ON` | Build GUI application |

### Conditional Compilation

| Define | Source | Description |
|--------|--------|-------------|
| `OPTIRAD_HAS_DCMTK` | `find_package(DCMTK QUIET)` | DICOM support available |
| `OPTIRAD_HAS_TBB` | `find_package(TBB QUIET)` | Parallel RT-STRUCT parsing |
| `OPTIRAD_DATA_DIR` | CMake cache variable | Path to `data/` directory |

### External Dependencies

| Dependency | Method | Required |
|------------|--------|----------|
| nlohmann/json v3.11.3 | `FetchContent` | Always |
| Dear ImGui | Vendored (`external/imgui/`) | GUI only |
| GoogleTest | `find_package` | Tests only |
| DCMTK | `find_package` (optional) | No |
| TBB | `find_package` (optional) | No |
| OpenMP | `find_package` (optional) | No |
| GLFW, GLEW, OpenGL, glm | `find_package` | GUI only |

## Threading Model

OptiRad uses three parallelism strategies:

1. **OpenMP** — Used in compute-heavy loops across most modules:
   - Dose calculation (per-beam, per-ray parallel iteration)
   - STF generation (parallel beam construction)
   - Voxel dilation (parallel neighbor search)
   - Dij matrix-vector products

2. **TBB** — Used for task-parallel RT-STRUCT parsing (when available):
   - `RTStructParser` parallelizes structure extraction across ROIs
    - `BodyContourGenerator` can parallelize per-slice BODY contour generation

3. **`std::thread`** — Used in the GUI for async long-running tasks:
    - Dose calculation, optimization, STF generation, and leaf sequencing run in background threads
   - Progress is communicated via `GuiAppState::taskProgress` and `taskStatus`
   - Cancellation via `std::atomic<bool> cancelFlag`

## Data Flow

### Central Data Container: `PatientData`

```
PatientData
├── Patient              → Demographics (name, ID)
├── Volume<int16_t>      → CT volume (Hounsfield Units)
├── Volume<double>        → Electron density (derived from CT)
└── StructureSet         → Collection of Structure objects
    └── Structure        → Name, type, color, contours, rasterized voxel indices
```

`PatientData` is created by `DicomImporter` and shared (via `std::shared_ptr`) across all pipeline stages.

### Dose Influence Matrix: `DoseInfluenceMatrix`

The Dij matrix maps bixel weights to voxel doses: $d = D_{ij} \cdot w$

- **Dimensions:** `numVoxels × numBixels` (typically millions × thousands)
- **Storage:** Sparse CSR after finalization
- **Key operations:**
  - `computeDose(weights)` → dose vector
  - `accumulateTransposeProduct(gradient, output)` → transpose for optimization gradient

### Optimization Loop

```
                    ┌─────────────────────┐
                    │   Initial weights    │
                    │   w = [1, 1, ..., 1] │
                    └─────────┬───────────┘
                              │
                    ┌─────────▼───────────┐
                    │  dose = Dij × w      │
                    └─────────┬───────────┘
                              │
                    ┌─────────▼───────────┐
                    │  Compute objectives  │
                    │  f(dose), ∇f(dose)   │
                    └─────────┬───────────┘
                              │
                    ┌─────────▼───────────┐
                    │  Gradient w.r.t. w   │
                    │  ∇w = Dij^T × ∇f    │
                    └─────────┬───────────┘
                              │
                    ┌─────────▼───────────┐
                    │  L-BFGS-B update     │
                    │  w ← w + α × d      │
                    │  (with bounds [0,∞)) │
                    └─────────┬───────────┘
                              │
                    ┌─────────▼───────────┐
                    │  Converged?          │──No──▶ Loop
                    └─────────┬───────────┘
                            Yes
                    ┌─────────▼───────────┐
                    │  OptimizationResult  │
                    └─────────────────────┘
```

## Module Documentation

For detailed documentation of each module, see:

- [Core Module](core.md) — Patient, Plan, Beam, Machine, Stf
- [Geometry Module](geometry.md) — Grid, Volume, Structure, coordinate systems
- [I/O Module](io.md) — DICOM import/export, machine loading
- [Dose Module](dose.md) — Dose engines, Dij matrix, ray tracing, plan analysis
- [Optimization Module](optimization.md) — L-BFGS-B optimizer, objective functions, constraints
- [Phase-Space Module](phsp.md) — IAEA format, beam sources
- [Segmentation Module](segmentation.md) — Automatic BODY contour generation from CT
- [Sequencing Module](sequencing.md) — Leaf sequencing and deliverable dose recomputation
- [Steering Module](steering.md) — STF generation
- [GUI Module](gui.md) — Application, panels, views, 3D renderers
- [Utilities Module](utils.md) — Logger, Config, Timer, math, interpolation
- [Acknowledgments](acknowledgments.md) — Upstream references and matRad attribution
