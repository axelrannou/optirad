# OptiRad

**Optimization of Radiotherapy Treatment Planning**

A modern, modular C++17 treatment planning system (TPS) for radiotherapy optimization with both CLI and GUI interfaces.

## Features

- **DICOM Import/Export**: Load CT images, RT Structure Sets, RT Plans, and RT Dose via DCMTK
- **Structure Management**: Import and visualize ROIs (targets, OARs, external) with contour rasterization
- **Dose Calculation**: Bortfeld SVD pencil beam engine with Siddon ray tracing, FFT lateral convolution, and electron density conversion
- **Dose Influence Matrix**: Sparse Dij (COO → CSR) with binary caching and deterministic filenames
- **Optimization**: L-BFGS-B optimizer with squared deviation/overdose/underdose and DVH-based objectives, NTO/hotspot control
- **Plan Analysis**: DVH computation, Dx%, VxGy, Conformity Index (CI), Homogeneity Index (HI), per-structure statistics
- **Phase-Space Beams**: IAEA phase-space file support (IAEA(NDS)-0484 format) for Monte Carlo beam sources
- **Steering (STF Generation)**: Target-aware ray generation with 3D margin expansion via morphological dilation
- **2D Visualization**: Interactive slice views (Axial/Sagittal/Coronal) with window/level, dose overlay, and contour rendering
- **3D Visualization**: GPU volume raycasting, marching cubes structure meshes, beam and phase-space particle rendering
- **GUI**: Modern Dear ImGui-based interface with 8 panels, 3 views, async task execution
- **CLI**: Interactive REPL with full planning pipeline (load → plan → STF → dose → optimize → analyze)
- **Machine Support**: Generic pencil beam machines (JSON) and Varian TrueBeam 6MV phase-space machines

## Quick Start

### Prerequisites

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake git \
    libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev libglew-dev libgtest-dev \
    libomp-dev libglm-dev

# Optional dependencies
sudo apt-get install -y libdcmtk-dev   # DICOM support
sudo apt-get install -y libtbb-dev     # Parallel RT-STRUCT parsing
```

### Build

```bash
git clone <repository-url>
cd optirad

# Setup GUI dependencies (Dear ImGui + OpenGL)
./scripts/setup_gui.sh

# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run tests
ctest --output-on-failure
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `OPTIRAD_BUILD_TESTS` | `ON` | Build GoogleTest unit tests |
| `OPTIRAD_BUILD_GUI` | `ON` | Build GUI application (requires OpenGL, GLFW, GLEW) |

## Usage

### Command-Line Interface (CLI)

The CLI provides an interactive REPL for the full treatment planning pipeline.

```bash
cd build

# Direct command execution
./apps/optirad_cli/optirad_cli load ../data/dicom/JOHN_DOE
./apps/optirad_cli/optirad_cli help

# Interactive mode (REPL)
./apps/optirad_cli/optirad_cli interactive
```

#### Available Commands

| Command | Description |
|---------|-------------|
| `load <dicom_dir>` | Load and inspect DICOM directory |
| `plan [options]` | Create a treatment plan |
| `generateStf` | Generate steering file (Generic machines) |
| `loadPhaseSpace [options]` | Load phase-space beam source (PhaseSpace machines) |
| `doseCalc [options]` | Calculate dose influence matrix (Dij) |
| `optimize [options]` | Run fluence optimization |
| `info` | Display current state information |
| `phsp-info` | Display phase-space source details |
| `interactive` | Enter interactive REPL mode |
| `help` | Show help message |

#### Plan Options

```
--mode <photons|protons>     Radiation mode (default: photons)
--machine <name>             Machine name (Generic, Varian_TrueBeam6MV)
--fractions <n>              Number of fractions (default: 30)
--gantry-start <deg>         Gantry start angle (default: 0)
--gantry-step <deg>          Gantry angle step (default: 4)
--gantry-stop <deg>          Gantry stop angle exclusive (default: 360)
--bixel-width <mm>           Bixel width (default: 7)
```

#### Dose Calculation Options

```
--dose-resolution <mm>       Dose grid resolution (default: 2.5)
--no-cache                   Disable Dij cache
```

#### Optimization Options

```
--max-iter <n>               Max iterations (default: 500)
--tolerance <val>            Convergence tolerance (default: 1e-5)
--target-dose <Gy>           Prescribed dose for targets (default: 60)
--oar-max-dose <Gy>          Max dose for OARs (default: 30)
```

#### Phase-Space Options

```
--collimator <deg>           Collimator angle (default: 0)
--couch <deg>                Couch angle (default: 0)
--max-particles <n>          Max particles per beam (default: 1000000)
--viz-sample <n>             Visualization sample per beam (default: 100000)
```

#### Example Workflow

```bash
./apps/optirad_cli/optirad_cli interactive

optirad> load ../data/dicom/JOHN_DOE
optirad> plan --machine Generic --gantry-step 40 --bixel-width 7
optirad> generateStf
optirad> doseCalc --dose-resolution 2.5
optirad> optimize --target-dose 60 --oar-max-dose 30 --max-iter 500
optirad> info
```

### Graphical User Interface (GUI)

```bash
cd build
./apps/optirad_gui/optirad_gui
```

#### GUI Panels

| Panel | Description |
|-------|-------------|
| **Patient Panel** | DICOM import dialog, file browser, patient demographics, structure list with visibility toggles |
| **Planning Panel** | Plan configuration (radiation mode, machine, gantry angles, bixel width, fractions), create plan / generate STF / calculate Dij |
| **STF Panel** | STF generation details, per-beam visibility toggles, progress tracking |
| **Phase-Space Panel** | Phase-space beam loading, per-beam visibility, particle statistics, energy histograms |
| **Optimization Panel** | Per-structure objective configuration (5 types), optimizer settings, NTO/hotspot control, async optimization |
| **Dose Stats Panel** | Plan analysis: stats table (min/max/mean, Dx%, VxGy, CI, HI), interactive DVH curves |
| **Beam Panel** | Beam information display |
| **Log Panel** | Log message viewer |

#### GUI Views

| View | Description |
|------|-------------|
| **Slice View** | 2D slice rendering (Axial/Sagittal/Coronal) with CT window/level, dose overlay (jet colormap), contour rendering |
| **3D View** | 3D viewport with orbit camera: GPU volume raycasting, structure meshes (marching cubes), beam rays, phase-space particles |
| **DVH View** | Dose-volume histogram curves per structure |

#### Window/Level Presets

| Preset | Window | Level | Use |
|--------|--------|-------|-----|
| Soft Tissue | 400 | 40 | General anatomy |
| Lung | 1500 | -600 | Lung parenchyma |
| Bone | 1800 | 400 | Skeletal structures |

## Data

### DICOM Support

| DICOM Type | SOP Class UID | Description |
|------------|---------------|-------------|
| CT Image Storage | 1.2.840.10008.5.1.4.1.1.2 | CT slices |
| RT Structure Set | 1.2.840.10008.5.1.4.1.1.481.3 | ROI contours |
| RT Plan | 1.2.840.10008.5.1.4.1.1.481.5 | Treatment plan |
| RT Dose | 1.2.840.10008.5.1.4.1.1.481.2 | Dose distribution |

### PatientData Structure

```
PatientData
├── Patient         // Demographics (name, ID)
├── CT Volume       // 3D array of HU values (int16_t)
├── ED Volume       // Electron density (double)
└── StructureSet    // Collection of ROIs
    └── Structures  // Individual contours + voxel masks
```

### Machine Configurations

Two machine types are supported:

- **Generic** — Pencil beam machine defined by a single JSON file (`data/machines/machine_photons_Generic.json`): SAD, SCD, energy, beta parameters, kernel data, fluence profile
- **Varian TrueBeam 6MV** — Phase-space machine defined by a folder (`data/machines/Varian_TrueBeam6MV/`): JSON config + multiple IAEA phase-space file pairs (`.IAEAheader` + `.IAEAphsp`)

### Dij Cache

Computed dose influence matrices are cached in `data/dij_cache/` with deterministic filenames based on patient name, beam count, bixel width, and resolution (e.g., `DOE^JOHN_90beams_bw7.0_res2.5mm.dij`).

### Sample Data

```
data/
├── dicom/
│   ├── JOHN_DOE/           # Sample DICOM dataset
│   ├── JANE_DOE/           # Sample DICOM dataset
│   └── JOHN_DOE_LUNG/      # Sample lung DICOM dataset
├── dij_cache/              # Cached Dij matrices (.dij binary)
├── machines/
│   ├── machine_photons_Generic.json
│   └── Varian_TrueBeam6MV/ # IAEA phase-space files
└── nrrd/                   # NRRD format data
```

## Configuration

### CT to Electron Density Conversion

Default piecewise linear HLUT (Hounsfield Lookup Table):

| HU Range | Electron Density | Tissue |
|----------|-----------------|--------|
| HU ≤ -1000 | ED = 0.0 | Air |
| -1000 < HU ≤ 0 | ED = 1.0 + HU/1000 | Lung / soft tissue |
| 0 < HU ≤ 100 | ED = 1.0 + 0.001 × HU | Soft tissue |
| HU > 100 | ED = 1.1 + 0.0005 × (HU - 100) | Bone |

## Project Structure

```
optirad/
├── CMakeLists.txt                      # Root CMake configuration
├── README.md
├── docs/                               # Technical documentation
│   ├── architecture.md                 # Global architecture overview
│   ├── core.md                         # Core module
│   ├── geometry.md                     # Geometry module
│   ├── io.md                           # I/O module
│   ├── dose.md                         # Dose calculation module
│   ├── optimization.md                 # Optimization module
│   ├── phsp.md                         # Phase-space module
│   ├── steering.md                     # Steering module
│   ├── gui.md                          # GUI module
│   └── utils.md                        # Utilities module
│
├── src/
│   ├── CMakeLists.txt
│   │
│   ├── core/                           # Core domain model
│   │   ├── Patient.hpp/cpp             # Patient demographics
│   │   ├── PatientData.hpp/cpp         # Central data container (CT + ED + structures)
│   │   ├── Plan.hpp/cpp                # Treatment plan
│   │   ├── Beam.hpp/cpp                # Beam definition (gantry/couch, rays, BEV/LPS)
│   │   ├── Ray.hpp                     # Ray/bixel definition
│   │   ├── Machine.hpp/cpp             # Linac parameters (Generic + PhaseSpace types)
│   │   └── Stf.hpp/cpp                 # Steering file (beam collection)
│   │
│   ├── geometry/                       # Geometry and spatial operations
│   │   ├── Grid.hpp/cpp                # 3D grid (voxel↔patient transforms)
│   │   ├── Volume.hpp/cpp              # Templated 3D volume (CTVolume, DoseVolume)
│   │   ├── Structure.hpp/cpp           # ROI with contours + rasterization
│   │   ├── StructureSet.hpp/cpp        # Collection of structures
│   │   ├── CoordinateSystem.hpp/cpp    # Patient↔world, voxel↔world transforms
│   │   ├── MathUtils.hpp               # Vec3, Mat3, rotation matrices, BEV↔LPS
│   │   └── VoxelDilation.hpp           # 3D morphological dilation (26-connectivity)
│   │
│   ├── io/                             # Input/Output
│   │   ├── IDataImporter.hpp           # Abstract importer interface
│   │   ├── IDataExporter.hpp           # Abstract exporter interface
│   │   ├── DicomImporter.hpp/cpp       # DICOM import (DCMTK)
│   │   ├── DicomExporter.hpp/cpp       # DICOM RT export
│   │   ├── RTStructParser.hpp/cpp      # RT-STRUCT parser (TBB parallel)
│   │   └── MachineLoader.hpp/cpp       # Machine JSON loader (auto-detect type)
│   │
│   ├── dose/                           # Dose calculation
│   │   ├── IDoseEngine.hpp             # Abstract dose engine interface
│   │   ├── DoseEngineFactory.hpp/cpp   # Factory: create engines by name
│   │   ├── DoseInfluenceMatrix.hpp/cpp # Sparse Dij matrix (COO → CSR)
│   │   ├── DoseMatrix.hpp/cpp          # 3D dose grid result
│   │   ├── DijSerializer.hpp/cpp       # Binary Dij serialization + caching
│   │   ├── SiddonRayTracer.hpp/cpp     # Siddon 1985 ray-voxel intersection
│   │   ├── SSDCalculator.hpp/cpp       # Source-to-surface distance
│   │   ├── RadDepthCalculator.hpp/cpp  # Radiological depth via ray tracing
│   │   ├── PlanAnalysis.hpp/cpp        # DVH, Dx%, VxGy, CI, HI statistics
│   │   ├── DoseCalcOptions.hpp         # Dose calculation options
│   │   ├── FFT2D.hpp                   # Cooley-Tukey 2D FFT (header-only)
│   │   ├── GridInterpolant2D.hpp       # Bilinear grid interpolation (header-only)
│   │   └── engines/
│   │       └── PencilBeamEngine.hpp/cpp  # Bortfeld SVD pencil beam
│   │
│   ├── optimization/                   # Optimization
│   │   ├── IOptimizer.hpp              # Abstract optimizer interface
│   │   ├── OptimizerFactory.hpp/cpp    # Factory: create optimizers by name
│   │   ├── ObjectiveFunction.hpp/cpp   # Base objective function
│   │   ├── Constraint.hpp/cpp          # Dose constraints (MinDose/MaxDose/MeanDose)
│   │   ├── objectives/
│   │   │   ├── SquaredDeviation.hpp/cpp   # (d - d_prescribed)²
│   │   │   ├── SquaredOverdose.hpp/cpp    # max(0, d - d_max)²
│   │   │   ├── SquaredUnderdose.hpp/cpp   # max(0, d_min - d)²
│   │   │   └── DVHObjective.hpp/cpp       # DVH-based (MIN_DVH / MAX_DVH)
│   │   └── optimizers/
│   │       └── LBFGSOptimizer.hpp/cpp  # L-BFGS-B with Wolfe line search
│   │
│   ├── phsp/                           # Phase-space beam sources
│   │   ├── IAEAHeaderParser.hpp/cpp    # IAEA header file parser
│   │   ├── IAEAPhspReader.hpp/cpp      # IAEA binary phase-space reader
│   │   ├── PhaseSpaceData.hpp/cpp      # Particle container with filtering/sampling
│   │   └── PhaseSpaceBeamSource.hpp/cpp # Full PSF pipeline (read → filter → transform)
│   │
│   ├── steering/                       # Steering file generation (header-only)
│   │   ├── StfProperties.hpp           # Lightweight STF properties
│   │   ├── IStfGenerator.hpp           # Abstract STF generator interface
│   │   └── PhotonIMRTStfGenerator.hpp  # Photon IMRT STF with target-aware rays
│   │
│   ├── gui/                            # Graphical User Interface
│   │   ├── Application.hpp/cpp         # Main GUI application
│   │   ├── Window.hpp/cpp              # GLFW window wrapper
│   │   ├── Renderer.hpp/cpp            # OpenGL frame management
│   │   ├── AppState.hpp                # Shared GUI application state
│   │   ├── panels/
│   │   │   ├── IPanel.hpp              # Abstract panel interface
│   │   │   ├── PatientPanel.hpp/cpp    # DICOM import, patient info, structures
│   │   │   ├── PlanningPanel.hpp/cpp   # Plan config, STF, Dij calculation
│   │   │   ├── StfPanel.hpp/cpp        # STF details, per-beam visibility
│   │   │   ├── PhaseSpacePanel.hpp/cpp # Phase-space loading, statistics
│   │   │   ├── OptimizationPanel.hpp/cpp # Objective config, optimizer settings
│   │   │   ├── DoseStatsPanel.hpp/cpp  # Plan analysis, DVH curves
│   │   │   ├── BeamPanel.hpp/cpp       # Beam information
│   │   │   └── LogPanel.hpp/cpp        # Log viewer
│   │   └── views/
│   │       ├── IView.hpp               # Abstract view interface
│   │       ├── SliceView.hpp/cpp       # 2D slice rendering with dose overlay
│   │       ├── View3D.hpp/cpp          # 3D viewport with orbit camera
│   │       ├── DVHView.hpp/cpp         # DVH curve plotting
│   │       └── renderers/
│   │           ├── VolumeRenderer.hpp/cpp       # GPU raycasting volume renderer
│   │           ├── StructureRenderer.hpp/cpp    # Marching cubes surface meshes
│   │           ├── BeamRenderer.hpp/cpp         # Beam ray visualization
│   │           ├── PhaseSpaceRenderer.hpp/cpp   # Phase-space particle rendering
│   │           ├── CubeRenderer.hpp/cpp         # Wireframe reference cube
│   │           ├── AxisLabels.hpp/cpp           # 3D axis labels
│   │           └── MarchingCubesTables.hpp      # Marching cubes lookup tables
│   │
│   └── utils/                          # Utilities
│       ├── Logger.hpp/cpp              # Thread-safe logging
│       ├── Config.hpp/cpp              # JSON configuration
│       ├── Timer.hpp/cpp               # High-resolution timer
│       ├── MathUtils.hpp/cpp           # Vec3 math operations
│       └── Interpolation.hpp/cpp       # Linear, bilinear, trilinear interpolation
│
├── apps/
│   ├── optirad_cli/
│   │   └── main.cpp                    # CLI with interactive REPL
│   └── optirad_gui/
│       └── main.cpp                    # GUI application entry point
│
├── tests/                              # GoogleTest unit tests (30 executables)
│   ├── core/                           # test_patient, test_plan, test_beam, test_ray
│   ├── geometry/                       # test_grid, test_volume, test_coordinate_system, test_indexing, test_voxel_dilation
│   ├── io/                             # test_dicom_importer, test_machine_loader
│   ├── dose/                           # test_dose_engine, test_dose_influence_matrix, test_dij_serializer, test_fft2d,
│   │                                   # test_grid_interpolant2d, test_siddon_ray_tracer, test_ssd_calculator,
│   │                                   # test_rad_depth_calculator, test_dij_matrad_comparison
│   ├── optimization/                   # test_lbfgs_optimizer
│   ├── phsp/                           # test_iaea_header_parser, test_iaea_phsp_reader, test_phase_space_data,
│   │                                   # test_phase_space_beam_source
│   ├── steering/                       # test_photon_imrt_stf_generator
│   └── utils/                          # test_logger, test_mathutils
│
├── external/                           # Third-party dependencies
│   └── imgui/                          # Dear ImGui (vendored)
│
├── scripts/
│   └── setup_gui.sh                    # Setup OpenGL + ImGui dependencies
│
└── data/                               # Runtime data
    ├── machines/                       # Machine definitions
    ├── dicom/                          # Sample DICOM datasets
    ├── dij_cache/                      # Cached Dij matrices
    └── nrrd/                           # NRRD format data
```

## Key Design Principles

### 1. Abstract Interfaces (Strategy Pattern)

- `IDoseEngine` — All dose engines implement this interface
- `IOptimizer` — All optimizers implement this interface
- `IDataImporter` / `IDataExporter` — Flexible I/O handling
- `IStfGenerator` — Steering file generation strategies
- `IPanel` / `IView` — Decoupled GUI components

### 2. Factory Pattern

- `DoseEngineFactory::create(name)` — Create dose engines by name
- `OptimizerFactory::create(name)` — Create optimizers by name

### 3. Pipeline Pattern

Full treatment planning pipeline:
`Load DICOM → Create Plan → Generate STF → Calculate Dij → Optimize → Analyze`

### 4. Modern C++17

- Smart pointers (`std::unique_ptr`, `std::shared_ptr`)
- `std::optional`, `std::variant` for type safety
- Templates for generic containers (`Volume<T>`)
- `constexpr`, structured bindings, `if constexpr`

### 5. Parallelism

- OpenMP for compute-heavy tasks (dose calculation, STF generation, dilation)
- TBB for parallel RT-STRUCT parsing
- `std::thread` for async GUI tasks with progress callbacks

## External Libraries

| Library | Purpose | Inclusion | Required |
|---------|---------|-----------|----------|
| [nlohmann/json](https://github.com/nlohmann/json) v3.11.3 | JSON parsing (machine files, config) | FetchContent | Yes |
| [Dear ImGui](https://github.com/ocornut/imgui) | Immediate mode GUI | Vendored (`external/imgui/`) | GUI only |
| [Google Test](https://github.com/google/googletest) | Unit testing framework | `find_package` | Tests only |
| [DCMTK](https://dicom.offis.de/dcmtk/) | DICOM import/export | `find_package` (optional) | No (`OPTIRAD_HAS_DCMTK`) |
| [TBB](https://github.com/oneapi-src/oneTBB) | Parallel RT-STRUCT parsing | `find_package` (optional) | No (`OPTIRAD_HAS_TBB`) |
| [OpenMP](https://www.openmp.org/) | Parallel computation | `find_package` (optional) | No |
| [GLFW](https://www.glfw.org/) | Window/input management | `find_package` | GUI only |
| [GLEW](http://glew.sourceforge.net/) | OpenGL extension loading | `find_package` | GUI only |
| [glm](https://github.com/g-truc/glm) | Math for 3D rendering | `find_package` | GUI only |
| OpenGL | Graphics rendering | System | GUI only |

## Documentation

See the [docs/](docs/) folder for detailed technical documentation:

- [Architecture Overview](docs/architecture.md) — System design, module dependencies, design patterns
- [Core Module](docs/core.md) — Patient, Plan, Beam, Machine, Stf
- [Geometry Module](docs/geometry.md) — Grid, Volume, Structure, coordinate systems
- [I/O Module](docs/io.md) — DICOM import/export, machine loading
- [Dose Module](docs/dose.md) — Dose engines, Dij matrix, ray tracing, plan analysis
- [Optimization Module](docs/optimization.md) — L-BFGS-B optimizer, objective functions, constraints
- [Phase-Space Module](docs/phsp.md) — IAEA format, beam sources
- [Steering Module](docs/steering.md) — STF generation
- [GUI Module](docs/gui.md) — Application, panels, views, 3D renderers
- [Utilities Module](docs/utils.md) — Logger, Config, Timer, math, interpolation
