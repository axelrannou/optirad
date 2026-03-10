# GUI Module (`optirad_gui`)

The GUI module provides a modern Dear ImGui-based graphical interface for interactive treatment planning with 2D/3D visualization, async task execution, and a panel-based layout.

**Library:** `optirad_gui`  
**Dependencies:** `optirad_core`, `optirad_io`, `optirad_geometry`, `optirad_optimization`, `optirad_phsp`, `imgui`, OpenGL, GLFW, GLEW, glm  
**Optional:** OpenMP

## Files

### Application Layer

| File | Description |
|------|-------------|
| `Application.hpp/cpp` | Main application class |
| `Window.hpp/cpp` | GLFW window wrapper |
| `Renderer.hpp/cpp` | OpenGL frame management |
| `AppState.hpp` | Shared application state |

### Panels (`panels/`)

| File | Description |
|------|-------------|
| `IPanel.hpp` | Abstract panel interface |
| `PatientPanel.hpp/cpp` | DICOM import, patient info, structure list |
| `PlanningPanel.hpp/cpp` | Plan configuration, STF, Dij calculation |
| `StfPanel.hpp/cpp` | STF details, per-beam visibility |
| `PhaseSpacePanel.hpp/cpp` | Phase-space beam loading and statistics |
| `OptimizationPanel.hpp/cpp` | Objective config, optimizer settings |
| `DoseStatsPanel.hpp/cpp` | Plan analysis, DVH curves |
| `BeamPanel.hpp/cpp` | Beam information display |
| `LogPanel.hpp/cpp` | Log message viewer |

### Views (`views/`)

| File | Description |
|------|-------------|
| `IView.hpp` | Abstract view interface |
| `SliceView.hpp/cpp` | 2D slice rendering with dose overlay |
| `View3D.hpp/cpp` | 3D viewport with orbit camera |
| `DVHView.hpp/cpp` | DVH curve plotting |

### 3D Renderers (`views/renderers/`)

| File | Description |
|------|-------------|
| `VolumeRenderer.hpp/cpp` | GPU raycasting volume renderer |
| `StructureRenderer.hpp/cpp` | Marching cubes surface meshes |
| `BeamRenderer.hpp/cpp` | Beam ray visualization |
| `PhaseSpaceRenderer.hpp/cpp` | Phase-space particle rendering |
| `CubeRenderer.hpp/cpp` | Wireframe reference cube |
| `AxisLabels.hpp/cpp` | 3D axis labels |
| `MarchingCubesTables.hpp` | Lookup tables for marching cubes |

## Architecture

### Application

Main entry point. Initializes GLFW window, ImGui context, and runs the main loop.

```cpp
class Application {
    bool init();       // Create window, init OpenGL, setup ImGui
    void run();        // Main loop: poll events → begin frame → render panels/views → end frame
    void shutdown();   // Cleanup
};
```

Internally holds:
- `GuiAppState` — Shared state container
- 3× `SliceView` — Axial, Sagittal, Coronal
- `View3D` — 3D viewport
- 8 Panels — Patient, Planning, STF, PhaseSpace, Optimization, DoseStats, Beam, Log

### GuiAppState

Central shared state container for the GUI. Extends the CLI's `AppState` with async task management and workflow queries.

```cpp
struct GuiAppState {
    // Pipeline data (shared with all panels/views)
    std::shared_ptr<PatientData>              patientData;
    std::shared_ptr<Plan>                     plan;
    std::shared_ptr<StfProperties>            stfProps;
    std::shared_ptr<Stf>                      stf;
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;
    std::shared_ptr<DoseInfluenceMatrix>      dij;
    std::shared_ptr<Grid>                     doseGrid;
    std::vector<double>                       optimizedWeights;
    std::shared_ptr<DoseMatrix>               doseResult;

    // Async task management
    std::atomic<bool> cancelFlag{false};
    std::string taskStatus;
    float taskProgress = 0.f;
    bool taskRunning = false;

    // Workflow queries
    bool dicomLoaded() const;
    bool planCreated() const;
    bool stfGenerated() const;
    bool phaseSpaceLoaded() const;
    bool dijComputed() const;
    bool optimizationDone() const;
    bool isPhaseSpaceMachine() const;

    // Cascade reset helpers
    void resetPlan();           // Resets plan + all downstream
    void resetStf();            // Resets STF + all downstream
    void resetPhaseSpace();
    void resetDij();
    void resetOptimization();
};
```

**Cascade resets:** When upstream data changes, all downstream data is invalidated. For example, `resetPlan()` also resets STF, Dij, and optimization results.

## Panels

### IPanel Interface

```cpp
class IPanel {
    virtual ~IPanel() = default;
    virtual void render(GuiAppState& state) = 0;
    virtual const char* getName() const = 0;
};
```

### PatientPanel

- File browser dialog for DICOM directory selection
- DICOM import trigger (`DicomImporter::importAll()`)
- Patient demographics display (name, ID)
- CT volume information (dimensions, spacing, origin, HU range)
- Structure list with:
  - Visibility checkboxes
  - Color indicators
  - Type labels (TARGET/OAR/EXTERNAL)
  - Contour/voxel counts

### PlanningPanel

- Radiation mode selector (photons)
- Machine selector (Generic, Varian_TrueBeam6MV)
- Gantry angle configuration (start, step, stop)
- Bixel width slider
- Number of fractions input
- **Create Plan** button → creates `Plan` with `StfProperties`
- **Generate STF** button → runs `PhotonIMRTStfGenerator::generateStf()` (async)
- **Calculate Dij** button → runs `PencilBeamEngine::calculateDij()` (async with progress)
- Dose resolution input
- Dij cache toggle

### StfPanel

- Beam count and total ray/bixel counts
- Per-beam summary table (gantry angle, ray count)
- Per-beam visibility toggles (for 3D view)
- Generation progress indicator

### PhaseSpacePanel

- Phase-space beam loading interface
- Collimator/couch angle inputs
- Max particles and visualization sample size
- Per-beam visibility toggles
- Particle count statistics (photons, electrons, positrons)
- Energy histogram display

### OptimizationPanel

- Per-structure objective configuration:
  - Objective type selector (SquaredDeviation, SquaredOverdose, SquaredUnderdose, DVH MIN, DVH MAX)
  - Dose/volume parameters
  - Weight sliders
- Optimizer settings:
  - Max iterations
  - Convergence tolerance
  - Max fluence
  - Memory size (L-BFGS)
- NTO/Hotspot control:
  - Prescription dose
  - Hotspot threshold
  - Hotspot penalty weight
- **Optimize** button → runs `LBFGSOptimizer::optimize()` (async)
- Progress bar with iteration count

### DoseStatsPanel

- Per-structure statistics table:
  - Min/Max/Mean/Std dose
  - D2, D5, D50, D95, D98
  - V20, V40, V50, V60
  - CI, HI (for targets)
- Interactive DVH curves (per-structure visibility toggles)

### BeamPanel

- Beam information display
- Summary of beam geometry

### LogPanel

- Scrolling log message viewer
- Severity-colored messages (info, warn, error, debug)

## Views

### SliceView

2D slice rendering with three orthogonal orientations.

**Features:**
- CT image display as OpenGL texture (grayscale)
- Window/level control (window width + center)
- Presets: Soft Tissue (W=400, L=40), Lung (W=1500, L=-600), Bone (W=1800, L=400)
- Slice navigation slider
- Structure contour overlay (colored lines)
- Dose overlay with jet colormap and adjustable opacity
- Coordinated crosshairs across views

**Orientations:**
| View | Plane | Axes |
|------|-------|------|
| Axial | X-Y | Left-Right × Anterior-Posterior |
| Sagittal | Y-Z | Anterior-Posterior × Superior-Inferior |
| Coronal | X-Z | Left-Right × Superior-Inferior |

### View3D

3D viewport with interactive orbit camera.

**Features:**
- Orbit camera (left-drag to rotate, scroll to zoom, middle-drag to pan)
- Multiple 3D renderers composited in a single viewport
- PIMPL idiom hides OpenGL implementation

**Uses renderers:**
- `VolumeRenderer` — CT volume raycasting
- `StructureRenderer` — Structure surface meshes
- `BeamRenderer` — Beam ray visualization
- `PhaseSpaceRenderer` — Phase-space particle scatter plot
- `CubeRenderer` — Reference wireframe cube
- `AxisLabels` — L/R, A/P, S/I orientation labels

### DVHView

Dose-volume histogram plot using ImGui drawing primitives.

- Per-structure DVH curves with structure colors
- Per-structure visibility toggles
- Dose axis (Gy) × Volume axis (%)

## 3D Renderers

### VolumeRenderer

GPU volume raycasting using OpenGL shaders.

**Algorithm:**
1. Upload CT volume as a 3D OpenGL texture (`GL_TEXTURE_3D`, `GL_R16I`)
2. Render a proxy cube (unit cube scaled to volume dimensions)
3. Fragment shader performs ray marching through the 3D texture
4. Apply transfer function (window/level) to map HU to opacity/color

### StructureRenderer

Structure surface mesh generation and rendering.

**Algorithm:**
1. For each visible structure, create a 3D binary mask from voxel indices
2. Run **marching cubes** algorithm on the mask to extract a triangulated isosurface
3. Upload vertex/index buffers to OpenGL
4. Render with structure color and basic lighting

The marching cubes lookup tables are in `MarchingCubesTables.hpp` (256 cases × 16 triangle indices).

### BeamRenderer

Beam geometry visualization.

- Source point → ray target lines (one per ray)
- Isocenter sphere marker
- Per-beam visibility control
- Color coding by beam index

### PhaseSpaceRenderer

Phase-space particle visualization as point clouds.

- Particle positions rendered as `GL_POINTS`
- Color mapping by energy or particle type
- Direction vectors (optional)
- Per-beam visibility control

### CubeRenderer

Wireframe reference cube showing the CT volume extent.

### AxisLabels

Text labels at the ends of coordinate axes (L/R, A/P, S/I) for patient orientation.

## Async Task Execution

Long-running computations (STF generation, dose calculation, optimization) run in background threads:

```cpp
// Example: async dose calculation
state.taskRunning = true;
state.taskStatus = "Calculating Dij...";
state.cancelFlag = false;

std::thread([&state]() {
    auto engine = DoseEngineFactory::create("PencilBeam");
    engine->setCancelFlag(&state.cancelFlag);
    engine->setProgressCallback([&](int beam, int total, const std::string& msg) {
        state.taskProgress = float(beam) / total;
        state.taskStatus = msg;
    });
    
    state.dij = std::make_shared<DoseInfluenceMatrix>(
        engine->calculateDij(*state.plan, *state.stf, *state.patientData, *state.doseGrid));
    
    state.taskRunning = false;
}).detach();
```

**Cancellation:** Panels check `state.taskRunning` to show progress and disable buttons. The cancel flag propagates to the engine via `IDoseEngine::setCancelFlag()`.

## Related Documentation

- [Architecture](architecture.md) — System overview and module dependencies
- [Core Module](core.md) — Data structures displayed and configured by panels
- [Dose Module](dose.md) — Dose calculation invoked from PlanningPanel
- [Optimization Module](optimization.md) — Optimization invoked from OptimizationPanel
- [Phase-Space Module](phsp.md) — Phase-space data rendered by PhaseSpaceRenderer
