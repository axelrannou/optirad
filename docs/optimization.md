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

## Constraints

### Constraint

Hard dose constraints for plan evaluation (not directly used in the optimization gradient, but checked for convergence/reporting).

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

    // NTO / Hotspot control
    void setPrescriptionDose(double dose);      // 0 = disabled
    void setHotspotThreshold(double threshold); // Fraction of Rx (default: 1.0)
    void setHotspotPenalty(double penalty);      // Penalty weight (default: 10000.0)
};
```

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
    
    // Optional: NTO hotspot penalty
    if hotspot control enabled:
        add penalty for voxels above threshold × Rx
    
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

#### NTO (Normal Tissue Objective) / Hotspot Control

When `prescriptionDose > 0`, the optimizer adds an additional penalty term for voxels exceeding the hotspot threshold:

$$f_\text{hotspot} = \lambda \sum_{v : d_v > \theta \cdot d_\text{Rx}} (d_v - \theta \cdot d_\text{Rx})^2$$

Where $\theta$ is `hotspotThreshold` (fraction of prescribed dose) and $\lambda$ is `hotspotPenalty`.

## Auto-Classification (CLI)

The CLI automatically classifies structures and creates objectives:

| Structure Type | Objectives Created |
|---------------|-------------------|
| **TARGET** | SquaredDeviation(Rx), MinDVH(D95 ≥ 0.95×Rx), MaxDVH(D2 ≤ 1.05×Rx) |
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
