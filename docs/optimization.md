# Optimization Module (`optirad_optimization`)

The optimization module provides fluence weight optimization for treatment planning: an abstract optimizer interface, L-BFGS-B implementation, objective functions, and dose constraints.

**Library:** `optirad_optimization`  
**Dependencies:** `optirad_dose`, `optirad_geometry`, `optirad_utils`  
**Optional:** OpenMP (parallel gradient computation)

## Files

| File | Description |
|------|-------------|
| `IOptimizer.hpp` | Abstract optimizer interface |
| `OptimizerFactory.hpp/cpp` | Factory for creating optimizers |
| `ObjectiveFunction.hpp/cpp` | Abstract base objective |
| `Constraint.hpp/cpp` | Dose constraints |
| `objectives/SquaredDeviation.hpp/cpp` | $(d - d_\text{Rx})^2$ |
| `objectives/SquaredOverdose.hpp/cpp` | $\max(0, d - d_\text{max})^2$ |
| `objectives/SquaredUnderdose.hpp/cpp` | $\max(0, d_\text{min} - d)^2$ |
| `objectives/DVHObjective.hpp/cpp` | DVH-based objective |
| `optimizers/LBFGSOptimizer.hpp/cpp` | L-BFGS-B optimizer |

## Interfaces

### IOptimizer

```cpp
struct OptimizationResult {
    std::vector<double> weights;
    double finalObjective;
    int iterations;
    bool converged;
};

class IOptimizer {
    virtual std::string getName() const = 0;

    virtual OptimizationResult optimize(
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        const std::vector<Constraint>& constraints) = 0;

    virtual void setMaxIterations(int maxIter) = 0;
    virtual void setTolerance(double tol) = 0;
};
```

### OptimizerFactory

```cpp
class OptimizerFactory {
    static std::unique_ptr<IOptimizer> create(const std::string& optimizerName);
};
```

Currently supports: `"LBFGS"`.

## Objective Functions

### ObjectiveFunction (Abstract Base)

All objective functions inherit from this base class, which manages weight, structure binding, and voxel index mapping.

```cpp
class ObjectiveFunction {
    virtual std::string getName() const = 0;
    virtual double compute(const std::vector<double>& dose) const = 0;
    virtual std::vector<double> gradient(const std::vector<double>& dose) const = 0;

    void setWeight(double weight);              // Importance weight (default: 1.0)
    double getWeight() const;
    void setStructure(const Structure* structure);
    void setVoxelIndices(const std::vector<size_t>& indices);   // Mapped dose grid indices

protected:
    const std::vector<size_t>& getActiveIndices() const;
    double m_weight = 1.0;
    const Structure* m_structure = nullptr;
    std::vector<size_t> m_mappedIndices;        // Indices into dose vector
};
```

**Voxel index mapping:** Structure voxel indices (on the CT grid) must be mapped to dose grid indices via `Grid::mapVoxelIndices()` before being set on the objective.

### SquaredDeviation

**Target objective.** Penalizes deviation from prescribed dose.

**Formula:** 

$$f = \frac{w}{N} \sum_{v \in \text{target}} (d_v - d_\text{Rx})^2$$

$$\nabla f_v = \frac{2w}{N} (d_v - d_\text{Rx})$$

```cpp
class SquaredDeviation : public ObjectiveFunction {
    void setPrescribedDose(double dose);        // Gy
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
};
```

### SquaredOverdose

**OAR objective.** Penalizes dose above a maximum threshold.

**Formula:**

$$f = \frac{w}{N} \sum_{v \in \text{OAR}} \max(0, d_v - d_\text{max})^2$$

$$\nabla f_v = \frac{2w}{N} \max(0, d_v - d_\text{max})$$

```cpp
class SquaredOverdose : public ObjectiveFunction {
    void setMaxDose(double dose);               // Gy
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
};
```

### SquaredUnderdose

**Target objective.** Penalizes dose below a minimum threshold.

**Formula:**

$$f = \frac{w}{N} \sum_{v \in \text{target}} \max(0, d_\text{min} - d_v)^2$$

$$\nabla f_v = \frac{-2w}{N} \max(0, d_\text{min} - d_v)$$

```cpp
class SquaredUnderdose : public ObjectiveFunction {
    void setMinDose(double dose);               // Gy
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
};
```

### DVHObjective

**DVH-based objective.** Penalizes based on dose-volume histogram constraints.

```cpp
enum class DVHObjective::Type { MIN_DVH, MAX_DVH };

class DVHObjective : public ObjectiveFunction {
    void setType(Type type);
    void setDoseThreshold(double dose);         // Gy
    void setVolumeFraction(double fraction);    // [0, 1]

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
};
```

| Type | Meaning | Use Case |
|------|---------|----------|
| `MIN_DVH` | At least `volumeFraction` of structure receives ≥ `doseThreshold` | D95 ≥ 57 Gy for targets |
| `MAX_DVH` | At most `volumeFraction` of structure receives ≥ `doseThreshold` | D2 ≤ 63 Gy (hotspot limit) |

Only the violating voxels count (the coldest / hottest voxels allowed by the volume are free):
$f = \frac{w}{N} \sum_{v \in \text{violating}} (d_v - d_\text{threshold})^2$.

**Priorities are volume-normalized.** Every objective, including the NTO, divides by its
structure voxel count $N$ (as in matRad and like Eclipse priorities), so a weight means the same
whatever the structure size, the objective type or the dose grid resolution. Weights of the order
of 100 for targets and the NTO are a good start; raise one objective in steps of ~50–100 if its
goal is not met.

## Constraints

### Constraint

Dose constraint value objects for plan evaluation. The current L-BFGS-B optimizer does **not**
enforce these objects or use them to determine convergence; it logs a warning if any are supplied.
Configure weighted dose objectives to penalize violations during optimization. Weight bounds and
DAO machine-deliverability limits are separate from these dose constraint objects.

```cpp
enum class ConstraintType { MinDose, MaxDose, MeanDose };

class Constraint {
    void setType(ConstraintType type);
    void setValue(double value);                 // Gy
    void setStructureName(const std::string& name);
    bool isSatisfied(double actualValue) const;
};
```

## L-BFGS-B Optimizer

### LBFGSOptimizer

Limited-memory BFGS with bound constraints.

**Inherits:** `IOptimizer`

```cpp
class LBFGSOptimizer : public IOptimizer {
    OptimizationResult optimize(
        const DoseInfluenceMatrix& dij,
        const std::vector<ObjectiveFunction*>& objectives,
        const std::vector<Constraint>& constraints) override;

    void setMaxIterations(int maxIter);         // default: 500
    void setTolerance(double tol);              // default: 1e-5

    // L-BFGS specific
    void setMemorySize(int m);                  // History size (default: 10)
    void setMaxFluence(double maxFluence);      // Upper bound (default: 10.0)
    void setVerbose(bool verbose);

};
```

#### GUI DAO progress

The **Optimization > DAO Progress** section displays live objective and total-MU
charts, resolution/iteration progress, maximum dose, and weighted penalty contribution bars.
Accepted-step counts and elapsed time are shown alongside expandable **Last sampled step**
details (gradient infinity norm, step size, and line-search iterations).
Choose a resolution from **History**, or leave **Latest resolution** selected to follow the run.
Each resolution has its own history because its dose model and penalty settings can change.
Charts use sampled-update indices, not uniformly spaced iterations. Log cost is enabled by
default (zero costs are floored at 1e-12 for plotting only).

Updates are published initially, after the first and every 10 accepted leaf/MU steps,
and after each resolution's feasibility repair. The final update is included even if no
steps were accepted. Histories remain attached to the resulting arc dose when switching
between saved dose entries; starting a new run clears the live display.

Bars are sorted by weighted cost and include dose objectives, the NTO, MU smoothing,
total MU, leaf order and leaf speed. Shares indicate cost contribution, **not clinical goal
attainment**. GUI snapshots add forward-dose/objective evaluations only when a progress
observer is attached. Terminal output retains its independent, quieter cadence. These
are optimization diagnostics, not live dose previews or final dose certification.

#### Terminal diagnostics

Both CLI and GUI optimization runs write diagnostics to standard output through `Logger`.
Launch the GUI from a terminal to see them. No additional flag is needed.

- **At startup:** each applied objective is listed by structure and readable type, with dose
  threshold in Gy, DVH volume percentage where applicable, weight, and active dose-grid voxel count.
  Objectives are explicitly identified as soft penalties, not hard clinical guarantees.
- **L-BFGS-B:** reports weight bounds, tolerance, smoothing and
  L1/L2 coefficients (including the derived relative L1 coefficient). Cost snapshots are emitted
  initially and finally; the compact iteration table continues at iteration 1 and every 10
  iterations. `setVerbose(false)` suppresses these snapshots and the iteration table, but retains
  configuration and warning logs.
- **DAO:** reports progressive-resolution stage transitions, seed creation or aperture/MU
  interpolation, seed MU calibration, actual dose-model backend, per-stage iteration budget,
  derived MU and mechanical penalty coefficients, and machine
  bounds. The initial cost breakdown is printed at the first stage and the final breakdown after
  feasibility repair at the last stage. Each stage prints a compact, aligned progress table on
  the first accepted step and every 50 accepted steps, with a header once per stage. Stages
  with no accepted steps do not print an empty table. Settings and summaries use short lines
  instead of a single long key/value list.
- **Stopping and delivery:** DAO reports the stop reason (budget, zero gradient, stationary
  bounds, line-search failure, or insufficient recent improvement), counts of history restarts
  and line-search retries, repair/snapping fallbacks, and final delivery MU, time, leaf speed and
  slowed control points.
  Final dose computation and pipeline completion are also logged.

DAO progress columns:

| Column | Meaning |
|--------|---------|
| `Iter` | Attempted iteration within this stage |
| `Acc` | Accepted steps within this stage |
| `Objective` | Current solver cost, including total-MU and mechanical penalties |
| `GradInf` | Pre-step infinity norm of the scaled-variable gradient (not a projected gradient) |
| `Step` | Accepted line-search step size; scientific notation keeps small steps visible |
| `LS` | Number of line-search trials for the accepted step |
| `MU` | Total MU of the accepted sequence |
| `Time(s)` | Wall time since the start of the stage's leaf refinement |

Each cost snapshot lists every objective's **already-weighted** cost and percentage of the total,
including the NTO, plus applicable regularization terms. DAO additionally separates
MU smoothing, total MU, leaf-order and leaf-speed penalties. Zero-cost terms remain visible.
Snapshot headers include the maximum dose on the computed grid.

**Interpretation:** cost shares describe optimization pressure, not clinical goal attainment,
gradient magnitude, or the percentage of voxels violating a goal. All objectives are
voxel-averaged, but a DVH objective only counts its violating voxels, so equal weights need
not produce equal costs. Use dose statistics/DVH for clinical evaluation.

DAO's `Objective` column includes total-MU and soft mechanical penalties. The summary's
`DoseCost` (`initialLeafObjective`/`finalLeafObjective`) includes dose objectives (including the NTO) and
MU smoothing only; it excludes total-MU and mechanical penalties. `Restarts` lists direction
restarts and line-search retries as `dir/LS`. Final snapshots describe
the repaired sequence, which may differ from the best unprojected iterate.

The two snapshots per optimization/refinement run each perform an additional forward-dose
evaluation and CPU evaluation of the objective costs, including on CUDA runs; progress lines
do not. No extra gradients or line-search trial logs are generated. CUDA solver totals and
snapshot totals may differ slightly due to floating-point reductions.

#### Algorithm

The L-BFGS-B optimizer uses the following approach:

1. **Initialization:** All weights set to 1.0
2. **Two-loop recursion:** Approximate inverse Hessian using last `m` gradient/step pairs
3. **Wolfe line search:** Find step size satisfying sufficient decrease and curvature conditions
4. **Bound projection:** Clamp weights to `[0, maxFluence]` after each update
5. **Convergence:** Stop when gradient norm < tolerance or max iterations reached

#### Optimization Loop

```
for iter = 1 to maxIterations:
    dose = Dij × weights
    
    // Compute total objective and gradient
    f = 0, ∇f = zeros(numVoxels)
    for each objective:
        f += objective.compute(dose)
        ∇f += objective.gradient(dose)
    
    // Transform gradient to weight space
    ∇w = Dij^T × ∇f
    
    // L-BFGS two-loop recursion → search direction
    d = -H_k × ∇w
    
    // Wolfe line search
    α = lineSearch(f, ∇w, d)
    
    // Update with bounds
    weights = clamp(weights + α × d, 0, maxFluence)
    
    if ||∇w|| < tolerance: converged = true; break
```

#### NTO (Normal Tissue Objective)

One Eclipse-style NTO (`NormalTissueObjective`, `NtoSettings`) controls the dose fall-off
outside the targets for both IMRT and VMAT/DAO. Each body voxel outside every target gets a
limit from its distance $x$ to the target border (Euclidean distance transform to the nearest
target voxel centre minus half a voxel):

$$L(x) = \begin{cases} s\,R_x & x \le x_0 \\ R_x\left(s\,e^{-k(x-x_0)} + e\,(1-e^{-k(x-x_0)})\right) & x > x_0 \end{cases}$$

with $x_0$ = distance from target border (default 1 mm), $s$ = start dose (105 %),
$e$ = end dose (60 %) and $k$ = fall-off (0.05 /mm, i.e. ~87 % of the way from start to end
10 mm beyond $x_0$). The cost is an upper objective averaged over the normal tissue voxels,
$\frac{\text{priority}}{N}\sum_v \max(0, d_v - L_v)^2$ (default priority 100), comparable with
target objectives of priority ~100. $R_x$ is the CLI `--target-dose`
(or plan prescription); in the GUI it is the *Prescription (Gy)* field (above the objectives) or, if 0, the highest
lower objective (MinDVH / Squared Deviation / Underdose) on a target — target upper objectives are
hot-spot limits and are not used. No prescription is ever assumed: `ObjectiveBuilder::build`,
`OptimizationPipeline` and `ArcPipeline` throw if it is not > 0, the CLI `optimize`/`arcSeq` stop
with an error when there is neither `--target-dose` nor a plan prescription, and the GUI disables
*Run* until one is known. The limit curve is logged at startup and previewed in the GUI.
Without a BODY/EXTERNAL structure the whole grid outside the targets is used (with a warning).

There is no separate global hot-spot penalty: as in Eclipse, control target hot spots with
Upper / Max DVH objectives on the target. A MaxDVH at volume 0% is a maximum dose (Eclipse Upper 0%);
at 100% it can never be violated and has no effect (likewise MinDVH at 0%). These cases, and
targets without any upper objective, are reported as warnings at startup and highlighted in the GUI. CLI: `--nto-priority`, `--nto-distance`,
`--nto-start`, `--nto-end`, `--nto-falloff`, `--no-nto` (for `optimize` and `arcSeq`).

The optimizer only pushes the dose up to the objective window it is given: a window such as
Lower 98% 62.7 Gy / Upper 0% 69.3 Gy leaves the target free anywhere between those doses. A
Squared Deviation to Rx makes the target homogeneous inside the window.

**Recommended start (GUI *Recommended* button, CLI target template).** Per target: Lower D98% ≥ Rx,
Upper 0% ≤ 105% Rx and Squared Deviation to Rx, priority 100 each; NTO priority 100 (1 mm, 105% →
60%, 0.05/mm); OAR protocol objectives; DAO MU smoothing 0.2. Rx is the GUI *Prescription (Gy)*
field, else the plan prescription (CLI: `--target-dose`, else the plan prescription); the GUI shows
the value and its source next to the button, and the template is neither built automatically nor
selectable until a prescription is set. On the lung test case (66 Gy, one
VMAT arc) this reaches PTV D95 ≈ 65 Gy, HI ≈ 0.045 and max ≈ 68.7 Gy (Eclipse plan: 64.8 Gy,
0.040, 68.3 Gy) with a lower body dose. Adjust from there: raise the priority of an objective that
is not met by ~50–100, or lower the NTO priority if the target coverage suffers.

## Auto-Classification (CLI)

The CLI automatically classifies structures and creates objectives:

| Structure Type | Objectives Created |
|---------------|-------------------|
| **TARGET** | MinDVH(D98 ≥ Rx), MaxDVH(D0 ≤ 1.05×Rx), SquaredDeviation(Rx), priority 100 |
| **OAR** | SquaredOverdose(maxDose) |
| **EXTERNAL** | (skipped) |

## Usage Example

```cpp
// Create optimizer
auto optimizer = OptimizerFactory::create("LBFGS");
optimizer->setMaxIterations(500);
optimizer->setTolerance(1e-5);

auto* lbfgs = dynamic_cast<LBFGSOptimizer*>(optimizer.get());
lbfgs->setMemorySize(10);
lbfgs->setMaxFluence(10.0);
lbfgs->setPrescriptionDose(60.0);

// Setup objectives
auto sqDev = std::make_unique<SquaredDeviation>();
sqDev->setPrescribedDose(60.0);
sqDev->setWeight(100.0);
sqDev->setStructure(targetStructure);
sqDev->setVoxelIndices(mappedTargetIndices);

auto sqOver = std::make_unique<SquaredOverdose>();
sqOver->setMaxDose(30.0);
sqOver->setWeight(50.0);
sqOver->setStructure(oarStructure);
sqOver->setVoxelIndices(mappedOARIndices);

std::vector<ObjectiveFunction*> objectives = {sqDev.get(), sqOver.get()};
std::vector<Constraint> constraints;

// Optimize
auto result = optimizer->optimize(dij, objectives, constraints);
// result.weights, result.finalObjective, result.iterations, result.converged
```

## Related Documentation

- [Dose Module](dose.md) — DoseInfluenceMatrix used for forward/gradient computation
- [Architecture](architecture.md) — Optimization loop diagram
