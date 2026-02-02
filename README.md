# optirad
OptiRad: Optimization of Radiotherapy Treatment Planning

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
│   │   ├── CMakeLists.txt
│   │   ├── Core.hpp               # Core module header
│   │   ├── Patient.hpp            # Patient data container
│   │   ├── Patient.cpp
│   │   ├── Plan.hpp               # Treatment plan
│   │   ├── Plan.cpp
│   │   ├── Beam.hpp               # Beam definition
│   │   ├── Beam.cpp
│   │   ├── Machine.hpp            # Machine/linac parameters
│   │   └── Machine.cpp
│   │
│   ├── io/                        # Input/Output module
│   │   ├── CMakeLists.txt
│   │   ├── IO.hpp                 # IO module header
│   │   ├── IDataImporter.hpp      # Abstract importer interface
│   │   ├── IDataExporter.hpp      # Abstract exporter interface
│   │   ├── DicomImporter.hpp      # DICOM import
│   │   ├── DicomImporter.cpp
│   │   ├── DicomExporter.hpp      # DICOM RT export
│   │   ├── DicomExporter.cpp
│   │   ├── NiftiImporter.hpp      # NIfTI support (optional)
│   │   └── NiftiImporter.cpp
│   │
│   ├── geometry/                  # Geometry and coordinate systems
│   │   ├── CMakeLists.txt
│   │   ├── Geometry.hpp           # Geometry module header
│   │   ├── CoordinateSystem.hpp   # World/patient coordinates
│   │   ├── CoordinateSystem.cpp
│   │   ├── Grid.hpp               # Dose grid definition
│   │   ├── Grid.cpp
│   │   ├── Volume.hpp             # 3D volume representation
│   │   ├── Volume.cpp
│   │   ├── Structure.hpp          # ROI/Structure definition
│   │   ├── Structure.cpp
│   │   ├── StructureSet.hpp       # Collection of structures
│   │   └── StructureSet.cpp
│   │
│   ├── dose/                      # Dose calculation engines
│   │   ├── CMakeLists.txt
│   │   ├── Dose.hpp               # Dose module header
│   │   ├── IDoseEngine.hpp        # Abstract dose engine interface
│   │   ├── DoseEngineFactory.hpp  # Factory for dose engines
│   │   ├── DoseEngineFactory.cpp
│   │   ├── DoseMatrix.hpp         # Dose matrix container
│   │   ├── DoseMatrix.cpp
│   │   ├── DoseInfluenceMatrix.hpp # Dij matrix for optimization
│   │   ├── DoseInfluenceMatrix.cpp
│   │   │
│   │   └── engines/               # Concrete dose engine implementations
│   │       ├── PencilBeamEngine.hpp
│   │       └── PencilBeamEngine.cpp
│   │   │
│   │   └── kernels/               # Dose kernels and data
│   │       ├── Kernel.hpp
│   │       ├── PhotonKernel.cpp
│   │       └── ProtonKernel.cpp
│   │
│   ├── optimization/              # Optimization module
│   │   ├── CMakeLists.txt
│   │   ├── Optimization.hpp       # Optimization module header
│   │   ├── IOptimizer.hpp         # Abstract optimizer interface
│   │   ├── OptimizerFactory.hpp   # Factory for optimizers
│   │   ├── OptimizerFactory.cpp
│   │   ├── ObjectiveFunction.hpp  # Base objective function
│   │   ├── ObjectiveFunction.cpp
│   │   ├── Constraint.hpp         # Constraint definitions
│   │   ├── Constraint.cpp
│   │   │
│   │   ├── objectives/            # Concrete objectives
│   │   │   ├── SquaredDeviation.hpp
│   │   │   ├── SquaredDeviation.cpp
│   │   │   ├── SquaredOverdose.hpp
│   │   │   ├── SquaredOverdose.cpp
│   │   │   ├── SquaredUnderdose.hpp
│   │   │   ├── SquaredUnderdose.cpp
│   │   │   ├── DVHObjective.hpp
│   │   │   ├── DVHObjective.cpp
│   │   │   ├── EUDObjective.hpp
│   │   │   └── EUDObjective.cpp
│   │   │
│   │   └── optimizers/            # Concrete optimizer implementations
│   │       ├── LBFGSOptimizer.hpp
│   │       └── LBFGSOptimizer.cpp
│   │
│   ├── gui/                       # Graphical User Interface (Dear ImGui)
│   │   ├── CMakeLists.txt
│   │   ├── GUI.hpp                # GUI module header
│   │   ├── Application.hpp        # Main application class
│   │   ├── Application.cpp
│   │   ├── Window.hpp             # GLFW/SDL window wrapper
│   │   ├── Window.cpp
│   │   ├── Renderer.hpp           # OpenGL/Vulkan renderer
│   │   ├── Renderer.cpp
│   │   │
│   │   ├── panels/                # ImGui panels/windows
│   │   │   ├── IPanel.hpp         # Abstract panel interface
│   │   │   ├── PatientPanel.hpp   # Patient info & structure list
│   │   │   ├── PatientPanel.cpp
│   │   │   ├── BeamPanel.hpp      # Beam configuration
│   │   │   ├── BeamPanel.cpp
│   │   │   ├── OptimizationPanel.hpp  # Objectives & constraints
│   │   │   ├── OptimizationPanel.cpp
│   │   │   ├── DoseStatsPanel.hpp # DVH, dose statistics
│   │   │   ├── DoseStatsPanel.cpp
│   │   │   ├── LogPanel.hpp       # Log/console output
│   │   │   └── LogPanel.cpp
│   │   │
│   │   ├── views/                 # 2D/3D visualization views
│   │   │   ├── IView.hpp          # Abstract view interface
│   │   │   ├── SliceView.hpp      # Axial/Sagittal/Coronal CT slices
│   │   │   ├── SliceView.cpp
│   │   │   ├── View3D.hpp         # 3D volume rendering / BEV
│   │   │   ├── View3D.cpp
│   │   │   ├── DVHView.hpp        # DVH plot (ImPlot)
│   │   │   ├── DVHView.cpp
│   │   │   ├── DoseProfileView.hpp # Dose profile curves
│   │   │   └── DoseProfileView.cpp
│   │   │
│   │   ├── widgets/               # Reusable ImGui widgets
│   │   │   ├── ColorMapWidget.hpp # Dose colormap selector
│   │   │   ├── SliderWidget.hpp   # Custom sliders
│   │   │   └── TransferFunctionWidget.hpp
│   │   │
│   │   └── rendering/             # OpenGL rendering utilities
│   │       ├── Shader.hpp         # Shader management
│   │       ├── Shader.cpp
│   │       ├── Texture.hpp        # CT/Dose textures
│   │       ├── Texture.cpp
│   │       ├── VolumeRenderer.hpp # 3D volume rendering
│   │       ├── VolumeRenderer.cpp
│   │       ├── SliceRenderer.hpp  # 2D slice rendering
│   │       └── SliceRenderer.cpp
│   │
│   └── utils/                     # Utilities and helpers
│       ├── CMakeLists.txt
│       ├── Utils.hpp              # Utils module header
│       ├── Logger.hpp             # Logging system
│       ├── Logger.cpp
│       ├── Config.hpp             # Configuration management
│       ├── Config.cpp
│       ├── Timer.hpp              # Performance timing
│       ├── Timer.cpp
│       ├── MathUtils.hpp          # Math helpers
│       ├── MathUtils.cpp
│       ├── Interpolation.hpp      # Interpolation functions
│       └── Interpolation.cpp
│
├── include/                       # Public headers (optional, for library use)
│   └── optirad/
│       └── OptiRad.hpp            # Main include header
│
├── apps/                          # Applications / executables
│   ├── CMakeLists.txt
│   ├── optirad_cli/               # CLI application
│   │   ├── CMakeLists.txt
│   │   └── main.cpp
│   └── optirad_gui/               # GUI application
│       ├── CMakeLists.txt
│       └── main.cpp
│
├── tests/                         # Unit and integration tests
│   ├── CMakeLists.txt
│   ├── core/
│   ├── io/
│   ├── geometry/
│   ├── dose/
│   └── optimization/
│
├── external/                      # Third-party dependencies (submodules/vendored)
│   ├── CMakeLists.txt
│   ├── eigen/                     # Linear algebra (header-only)
│   ├── spdlog/                    # Logging
│   ├── nlohmann_json/             # JSON parsing
│   └── dcmtk/                     # DICOM toolkit
│
├── data/                          # Sample data and resources
│   ├── machines/                  # Machine definition files
│   ├── kernels/                   # Dose kernels
│   └── test_data/                 # Test datasets
│
└── scripts/                       # Build and utility scripts
    ├── build.sh
    └── run_tests.sh
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
- Concepts (C++20) for template constraints

## Recommended External Libraries

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
| glad/glew | OpenGL loader | MIT |

## Getting Started

```bash
# Clone and setup
git clone <repository>
cd optirad
mkdir build && cd build

# Configure and build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)

# Run tests
ctest --output-on-failure
