# GUI Module (`optirad_gui`)

The GUI module provides a modern Dear ImGui-based graphical interface for interactive treatment planning with 2D/3D visualization, async task execution, and a panel-based layout.

**Library:** `optirad_gui`  
**Dependencies:** `optirad_workflow`, `optirad_core`, `optirad_io`, `optirad_geometry`, `optirad_optimization`, `optirad_phsp`, `imgui`, OpenGL, GLFW, GLEW, glm  
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
| `DoseStatsPanel.hpp/cpp` | Plan analysis statistics table |
| `DVHPanel.hpp/cpp` | Dedicated DVH panel with comparison overlay |
| `LeafSequencingPanel.hpp/cpp` | Leaf sequencing controls and deliverable-dose results |
| `BeamPanel.hpp/cpp` | Beam information display |
| `LogPanel.hpp/cpp` | Log message viewer |

### Views (`views/`)

| File | Description |
|------|-------------|
| `IView.hpp` | Abstract view interface |
| `SliceView.hpp/cpp` | 2D slice rendering with dose overlay |
| `View3D.hpp/cpp` | 3D viewport with orbit camera |
| `DVHView.hpp/cpp` | DVH curve plotting |
| `BevView.hpp/cpp` | Beam's Eye View with fluence and aperture visualization |

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
- `DVHView` and `BevView` — Specialized dose-analysis and beam-eye displays
- 10 Panels — Patient, Planning, STF, PhaseSpace, Optimization, DoseStats, DVH, LeafSequencing, Beam, Log

### GuiAppState

Central shared state container for the GUI. It builds on the shared workflow-state model and adds async task handling, multi-dose bookkeeping, and GUI refresh state.

```cpp
struct GuiAppState {
  // Core workflow state
  std::shared_ptr<PatientData> patientData;
  std::shared_ptr<Plan> plan;
  std::shared_ptr<StfProperties> stfProps;
  std::shared_ptr<Stf> stf;
  std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;
  std::shared_ptr<DoseInfluenceMatrix> dij;
  std::shared_ptr<Grid> doseGrid;
  std::shared_ptr<Grid> computeGrid;

  // Optimization and sequencing state
  std::vector<double> optimizedWeights;
  std::shared_ptr<DoseMatrix> doseResult;
  std::vector<LeafSequenceResult> leafSequences;
  std::vector<double> deliverableWeights;
  std::shared_ptr<DoseMatrix> deliverableDoseResult;
  std::vector<StructureDoseStats> deliverableStats;

  // Multi-dose management and caches
  DoseManager doseManager;
  std::unordered_map<int, std::vector<double>> optWeightsCache;
  std::unordered_map<int, LeafSeqCacheEntry> seqCache;

  // Async task management
  std::atomic<bool> cancelFlag{false};
  std::string taskStatus;
  float taskProgress = 0.f;
  bool taskRunning = false;
  std::atomic<bool> optimizationJustFinished{false};

  // Workflow queries and reset helpers
  bool dicomLoaded() const;
  bool planCreated() const;
  bool stfGenerated() const;
  bool phaseSpaceLoaded() const;
  bool dijComputed() const;
  bool optimizationDone() const;
  bool leafSequenceDone() const;
  bool doseAvailable() const;
  bool isPhaseSpaceMachine() const;
  void syncSelectedDose();
  void resetPlan();
  void resetStf();
  void resetPhaseSpace();
  void resetDij();
  void resetLeafSequence();
  void resetOptimization();
  void resetAllDoses();
};
```

**Cascade resets:** When upstream data changes, downstream pipeline state is invalidated while the selected display dose can remain available through `DoseManager`.

**VMAT (DAO only):** a VMAT plan is always optimized by direct aperture optimization. The Planning panel's Dose Calculation section is replaced by a few dose settings (grid, backend, threads); no Dij, fluence optimization or sequencing step is needed. The Optimization panel holds the objectives, the NTO, the DAO settings, the *Run DAO* button and the live DAO progress; the Leaf Sequencing panel only shows the resulting apertures and control points.

**Multi-dose support:** The GUI can keep multiple dose maps alive at once, including imported dose, optimization results, and deliverable dose after leaf sequencing. Panels use `DoseManager` versioning and cached per-dose statistics so switching views does not force recomputation every time.

## Panels

### IPanel Interface

```cpp
class IPanel {
    virtual ~IPanel() = default;
  virtual std::string getName() const = 0;
  virtual void update() = 0;
  virtual void render() = 0;

  void setVisible(bool visible);
  bool isVisible() const;
};
```

### PatientPanel

- File browser dialog for DICOM directory selection
- DICOM import trigger (`DicomImporter::importAll()`)
- Imported RT Dose is selected and available to DVH and Dose Statistics immediately
  when CT and structures are present; no plan or optimization is required.
- A successful patient import replaces the previous patient's plan, doses, and
  analysis caches. A failed import preserves the current patient.
- Patient demographics display (name, ID)
- CT volume information (dimensions, spacing, origin, HU range)
- Structure list with:
  - Visibility checkboxes
  - Color indicators
  - Type labels (TARGET/OAR/EXTERNAL)
  - Contour/voxel counts

### PlanningPanel

- Radiation mode selector (photons)
- Machine selector (Generic, HD120, Varian_TrueBeam6MV; default HD120)
- Delivery defaults to VMAT with a clockwise 0 to 360 degree arc, 10 degree sectors,
  2 degree control point spacing, and 3 mm bixels (also the CLI planning defaults)
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
  - Weight sliders (volume-normalized priorities: ~100 for targets and the NTO is a good start)
- **Prescription (Gy)** field, with the template Rx and its source shown next to it ("Prescription
  field" or "plan prescription"); 0 = use the plan prescription
- **Recommended** button: replaces the objectives and NTO by the tested template (per PTV Lower
  D98% ≥ Rx, Upper 0% ≤ 105% Rx, Squared Deviation Rx, priority 100; OAR protocol; NTO 100).
  Disabled, and the template is not loaded for a new patient, until a prescription is known
  (no dose is assumed); *Run* is disabled without a prescription
- Default OAR objectives (also restored by **Recommended**, matched by exact structure name):

  | Structure | Objective | Dose (Gy) | Volume (%) | Weight |
  |-----------|-----------|-----------|------------|--------|
  | Canal_Medullaire | MaxDVH | 45 | 0 | 150 |
  | Poumon_D, Poumon_G (each) | MaxDVH | 20 | 30 | 150 |
  | Coeur | MaxDVH | 40 | 30 | 50 |
  | Oesophage | MaxDVH | 50 | 30 | 50 |

  These are soft optimization penalties, not guaranteed hard constraints. Lung mean dose
  <= 20 Gy remains a separate plan-evaluation criterion; the V20 objective does not enforce it.
  Evaluate the requested strict V20/V40/V50 < 30% criteria after optimization: MaxDVH targets
  an allowed 30% volume, rather than a strict inequality. Per-lung objectives do not replace
  evaluating a combined-lung volume when required by the protocol. Cord and lung limits take
  precedence over heart/oesophagus goals during plan review; weights alone do not guarantee
  this hierarchy. These curated GUI OAR presets do not change the CLI's objective selection.
- Optimizer settings:
  - Max iterations
  - Convergence tolerance
  - Max fluence
  - Memory size (L-BFGS)
- Normal Tissue Objective (Eclipse-style, one for IMRT and VMAT): Use NTO, the Rx it uses (Prescription
  field, else the highest lower target objective), priority, distance from target border (mm, default 1),
  start dose (% Rx), end dose (% Rx), fall-off (/mm), and a preview of the limit curve (0–50 mm)
- VMAT plans: Direct Aperture Optimization settings (iterations per resolution, first control
  point spacing, MU penalty, relative MU smoothing (default 0.2), integrate leaf motion), **Run DAO** button, status and
  the live DAO Progress charts/bars (`DaoOptimizationView`)
- IMRT plans: **Optimize** button → runs `LBFGSOptimizer::optimize()` (async), progress bar

### DoseStatsPanel

- Per-structure statistics table:
  - Min/Max/Mean/Std dose
  - D2, D5, D50, D95, D98
  - V20, V40, V50, V60
  - CI, HI (for targets)
  - Included CT voxel count, total structure voxel count, and coverage percentage

The Dose Statistics panel, DVH panel, and both CSV exports use the **same CT
structure voxel centres inside every loaded dose grid**. Dose values are
trilinearly interpolated at these centres; native dose-grid resolutions may
differ. The existing half-voxel boundary support is retained. Outside-grid dose
is unknown, not zero-filled. Genuine zero-dose samples inside the common region
remain included.

Coverage below 100% means all reported metrics describe only the common part of
the organ, not the whole organ. Structures with no common samples are omitted
with a log warning. Adding or removing a dose map invalidates the statistics
cache and changes the common region for every loaded map, even if it is not
selected for the comparison overlay.

CSV `NumVoxels` is the included CT sample count. `TotalStructureVoxels` and
`Coverage(%)` are appended to report the full mask size and included percentage.
Existing CSVs must be regenerated from the dose maps; their aggregate statistics
cannot be corrected by editing voxel counts. Standalone pipeline/CLI analysis
retains its single-dose coverage unless comparison grids are supplied explicitly.

### DVHPanel

- Displays DVH curves for the selected dose map
- Supports an optional comparison dose overlay
- Tracks dose-version changes and recomputes only when needed
- Uses per-structure visibility toggles and configurable maximum dose range

### LeafSequencingPanel

- IMRT: activates after optimization completes, configures intensity levels and minimum segment
  MU, runs `LeafSequencingPipeline` asynchronously and displays per-beam segments
- VMAT: results only (summary, DAO accepted steps and objective, speeds, control point table);
  DAO is configured and run from the Optimization panel
- Integrates with `DoseManager` so deliverable dose becomes selectable like any other dose map

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

### BevView

Beam's Eye View for inspection of sequenced and optimized beam fluence.

- Displays a fluence heatmap for the selected beam
- Overlays MLC leaf positions for a selected segment or the combined fluence
- Optionally draws projected patient structure contours in BEV coordinates
- Supports zoom and pan similar to `SliceView`
- Refreshes automatically when dose or sequencing selection changes

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

Long-running computations (STF generation, dose calculation, optimization, leaf sequencing) run in background threads:

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
- [Sequencing Module](sequencing.md) — Leaf sequencing data shown in `LeafSequencingPanel` and `BevView`
