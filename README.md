# optirad
OptiRad: Optimization of Radiotherapy Treatment Planning

A modern, modular treatment planning system (TPS) for radiotherapy optimization with both CLI and GUI interfaces.

## Features

- **DICOM Import**: Load CT images, RT structures, RT plans, and RT dose
- **Structure Management**: Import and visualize ROIs (targets, OARs)
- **Dose Calculation**: Pencil beam dose engine with electron density conversion
- **Optimization**: L-BFGS-B optimizer with multiple objective functions
- **Visualization**: Interactive 2D slice views (Axial/Sagittal/Coronal) with window/level control
- **GUI**: Modern ImGui-based interface with real-time interaction
- **CLI**: Command-line tools for batch processing and scripting

## Quick Start

### Prerequisites

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake git \
    libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev \
    libomp-dev

# Optional: DCMTK for DICOM support
sudo apt-get install -y libdcmtk-dev
```

### Build

```bash
# Clone repository
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

## Usage

### Command-Line Interface (CLI)

The CLI application provides a simple interface for batch processing and inspection.

#### Load and Inspect DICOM Data

```bash
# Load DICOM directory
./apps/optirad_cli/optirad_cli load /path/to/dicom/folder
```

#### Available Commands

```bash
# Show help
./optirad_cli help

# Load DICOM directory
./optirad_cli load <dicom_directory>
```

### Graphical User Interface (GUI)

The GUI provides an interactive visualization and planning environment.

#### Launch GUI

```bash
./apps/optirad_gui/optirad_gui
```

#### GUI Features

**1. Import DICOM Data**
- Click **"Import DICOM Directory"** button in the Patient Panel
- Browse to your DICOM folder
- The system will automatically detect and load:
  - CT image series (sorted by position)
  - RT Structure Set with contours
  - RT Plan (if available)
  - RT Dose (if available)

**2. View Patient Information**
- **Patient Panel** (left side): Shows patient demographics, CT dimensions, and structure list
- Toggle structure visibility with checkboxes
- Color indicators show structure colors

**3. Slice Visualization**
Three orthogonal views are available:
- **Axial View**: X-Y plane (looking down through patient)
- **Sagittal View**: Y-Z plane (side view)
- **Coronal View**: X-Z plane (front view)

**4. Window/Level Control**
Each slice view includes:
- **Slice Slider**: Navigate through slices
- **Window Width**: Adjust contrast (HU range)
- **Window Center**: Adjust brightness (center HU)
- **Presets**: Quick settings for different tissues:
  - **Soft Tissue**: W=400, L=40 (general anatomy)
  - **Lung**: W=1500, L=-600 (lung parenchyma)
  - **Bone**: W=1800, L=400 (skeletal structures)

## Data Format

### DICOM Support

OptiRad supports standard DICOM-RT objects:

| DICOM Type | SOP Class UID | Description |
|------------|---------------|-------------|
| CT Image Storage | 1.2.840.10008.5.1.4.1.1.2 | CT slices |
| RT Structure Set Storage | 1.2.840.10008.5.1.4.1.1.481.3 | ROI contours |
| RT Plan Storage | 1.2.840.10008.5.1.4.1.1.481.5 | Treatment plan |
| RT Dose Storage | 1.2.840.10008.5.1.4.1.1.481.2 | Dose distribution |

### PatientData Structure

Internally, data is organized similar to matRad's `ct` + `cst` structure:

```cpp
PatientData
├── Patient         // Demographics (name, ID)
├── CT Volume       // 3D array of HU values (int16_t)
├── ED Volume       // Electron density (double)
└── StructureSet    // Collection of ROIs
    └── Structures  // Individual contours + voxel masks
```

### Example Data Layout

```
data/dicom/PATIENT_NAME/
├── CT_SERIES_UID/
│   ├── slice_001.dcm
│   ├── slice_002.dcm
│   └── ...
├── RS_SERIES_UID/
│   └── rtstruct.dcm
├── RTPLAN_SERIES_UID/
│   └── rtplan.dcm
└── RTDOSE_SERIES_UID/
    └── rtdose.dcm
```

## Configuration

### CT to Electron Density Conversion

Default piecewise linear HLUT (Hounsfield Lookup Table):

```cpp
HU ≤ -1000  → ED = 0.0     (air)
-1000 < HU ≤ 0   → ED = 1.0 + HU/1000  (lung/soft tissue)
0 < HU ≤ 100     → ED = 1.0 + 0.001*HU (soft tissue)
HU > 100         → ED = 1.1 + 0.0005*(HU-100) (bone)
```

## Project Structure / Architecture

```
optirad/
├── CMakeLists.txt                 # Root CMake configuration
├── README.md
├── docs/                          # Documentation
│   └── architecture.md
│
├── src/
│   ├── CMakeLists.txt
│   │
│   ├── core/                      # Core abstractions and base classes
│   │   ├── Patient.hpp            # Patient demographics
│   │   ├── PatientData.hpp        # Main data container (ct + cst)
│   │   ├── Plan.hpp               # Treatment plan
│   │   ├── Beam.hpp               # Beam definition
│   │   └── Machine.hpp            # Linac parameters
│   │
│   ├── io/                        # Input/Output module
│   │   ├── DicomImporter.hpp      # DICOM import with DCMTK
│   │   ├── DicomExporter.hpp      # DICOM RT export
│   │   └── IDataImporter.hpp      # Abstract interface
│   │
│   ├── geometry/                  # Geometry and coordinate systems
│   │   ├── Grid.hpp               # 3D dose grid
│   │   ├── Volume.hpp             # Templated 3D volume
│   │   ├── Structure.hpp          # ROI with contours
│   │   └── StructureSet.hpp       # Collection of structures
│   │
│   ├── dose/                      # Dose calculation engines
│   │   ├── IDoseEngine.hpp        # Abstract dose engine
│   │   ├── DoseInfluenceMatrix.hpp # Dij matrix (sparse)
│   │   └── engines/               # Concrete implementations
│   │       └── PencilBeamEngine.hpp
│   │
│   ├── optimization/              # Optimization module
│   │   ├── IOptimizer.hpp         # Abstract optimizer
│   │   ├── ObjectiveFunction.hpp  # Base objective
│   │   ├── objectives/            # Concrete objectives
│   │   │   ├── SquaredDeviation.hpp
│   │   │   ├── SquaredOverdose.hpp
│   │   │   └── SquaredUnderdose.hpp
│   │   └── optimizers/            # Concrete optimizers
│   │       └── LBFGSOptimizer.hpp # L-BFGS-B implementation
│   │
│   ├── gui/                       # Graphical User Interface
│   │   ├── Application.hpp        # Main GUI application
│   │   ├── panels/                # UI panels
│   │   │   ├── PatientPanel.hpp   # Patient info & structures
│   │   │   └── OptimizationPanel.hpp
│   │   └── views/                 # Visualization views
│   │       ├── SliceView.hpp      # 2D slice rendering
│   │       └── DVHView.hpp        # DVH plots
│   │
│   └── utils/                     # Utilities
│       ├── Logger.hpp             # Logging system
│       └── MathUtils.hpp          # Math helpers
│
├── apps/                          # Applications
│   ├── optirad_cli/               # CLI application
│   │   └── main.cpp
│   └── optirad_gui/               # GUI application
│       └── main.cpp
│
├── tests/                         # Unit tests
│   ├── io/test_dicom_importer.cpp
│   └── optimization/test_lbfgs_optimizer.cpp
│
├── external/                      # Third-party dependencies
│   └── imgui/                     # Dear ImGui (auto-downloaded)
│
└── scripts/                       # Build scripts
    └── setup_gui.sh               # Setup OpenGL + ImGui
```

## Key Design Principles

### 1. **Abstract Interfaces (Strategy Pattern)**
- `IDoseEngine` - All dose engines implement this interface
- `IOptimizer` - All optimizers implement this interface
- `IDataImporter/IDataExporter` - Flexible IO handling
- `IPanel/IView` - Decoupled GUI panels and views

### 2. **Factory Pattern**
- `DoseEngineFactory` - Create dose engines by name/type
- `OptimizerFactory` - Create optimizers by name/type

### 3. **Dependency Injection**
- Components receive dependencies via constructor
- Enables easy testing and swapping implementations

### 4. **Modern C++ (C++17/20)**
- Smart pointers (`std::unique_ptr`, `std::shared_ptr`)
- `std::optional`, `std::variant` for type safety
- Templates for generic containers (`Volume<T>`)

## External Libraries

| Library | Purpose | License |
|---------|---------|---------|
| Eigen | Linear algebra, matrices | MPL2 |
| DCMTK | DICOM import/export | BSD |
| spdlog | Fast logging | MIT |
| nlohmann/json | JSON config files | MIT |
| Google Test | Unit testing | BSD |
| Dear ImGui | Immediate mode GUI | MIT |
| ImPlot | Plotting for ImGui | MIT |
| GLFW | Window/input management | Zlib |
| OpenGL | Graphics rendering | - |
