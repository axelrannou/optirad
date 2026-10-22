# Sequencing Module (`optirad_sequencing`)

The sequencing module converts optimized beamlet fluence into step-and-shoot MLC apertures and computes the associated deliverable dose.

**Library:** `optirad_sequencing`
**Dependencies:** `optirad_core`, `optirad_dose`

## Files

| File | Description |
|------|-------------|
| `LeafSequencer.hpp/cpp` | Step-and-shoot leaf sequencing from fluence maps |
| `DeliverableDoseCalculator.hpp/cpp` | Reconstruct deliverable weights and dose from sequenced apertures |

## Main Types

### Aperture and Sequencing Result

The sequencing module relies on the core types declared in `core/Aperture.hpp`.

- `Aperture` stores leaf bank positions and segment MU.
- `LeafSequenceResult` stores per-beam sequencing output, fidelity metrics, quantized fluence, and optional leaf-pair geometry.
- `LeafSequencerOptions` controls quantization levels, minimum segment MU, and leaf position snapping resolution.

## Classes

### LeafSequencer

Converts a 2D `FluenceMap` into a set of step-and-shoot apertures.

```cpp
class LeafSequencer {
public:
    static LeafSequenceResult sequenceBeam(
        const FluenceMap& fluence,
        const MachineGeometry& mlc,
        const LeafSequencerOptions& opts);
};
```

#### Sequencing Workflow

The implementation performs the following high-level steps:

1. Build leaf-pair boundaries from machine MLC geometry.
2. Determine which leaf pairs overlap the fluence extent.
3. Resample bixel-grid fluence onto the physical leaf-pair grid.
4. Quantize the fluence into a fixed number of intensity levels.
5. Back-project the quantized fluence to the bixel grid for later dose reconstruction.
6. Extract apertures level by level.
7. Snap leaf positions to the configured mechanical resolution.
8. Merge or remove segments below the minimum MU threshold when requested.
9. Compute sequencing quality metrics such as total MU and fluence fidelity.

The result stores both the deliverable aperture sequence and the quantized fluence representation needed for later dose evaluation and BEV visualization.

### DeliverableDoseCalculator

Computes effective bixel weights and dose from the sequenced representation.

```cpp
class DeliverableDoseCalculator {
public:
    static std::vector<double> computeDeliverableWeights(
        const std::vector<LeafSequenceResult>& sequences,
        const Stf& stf);

    static std::shared_ptr<DoseMatrix> computeDeliverableDose(
        const DoseInfluenceMatrix& dij,
        const std::vector<double>& deliverableWeights,
        const Grid& doseGrid);
};
```

This allows comparison between:

- the ideal optimized fluence,
- the sequenced, machine-deliverable fluence,
- and the resulting deliverable dose distribution.

## Workflow Integration

The sequencing module is typically used after optimization:

```text
optimized weights -> FluenceMap -> LeafSequencer -> deliverable weights -> dose recomputation
```

That workflow is orchestrated in the shared pipeline layer by `core/workflow/LeafSequencingPipeline`.

## Notes

- The module currently targets step-and-shoot IMRT delivery.
- Mixed-width MLCs are supported through the machine geometry leaf-width description.
- Sequencing quality is reported using fluence fidelity metrics to help compare the deliverable result against the optimized fluence.
- Several sequencing conventions were chosen to align comparisons with matRad-style workflows. See `docs/acknowledgments.md`.

## VMAT arc sequencing

Delivery mode is chosen at planning time (`SMLC` or `VMAT`, `Plan::setDeliveryMode`). For VMAT the arc is defined by start angle, stop angle, direction and couch angle (`ArcSpec`), plus an optimization-beam spacing (default 20 deg) and a control-point spacing (default 4 deg) (`VmatSettings`).

Workflow:
1. Fluence optimization runs on one beam per arc sector (sector centres) as usual.
2. `ArcSequencer::sequence` turns each sector fluence map into one aperture per control point.
   - **Sweep** (default; needs an interdigitating MLC and >= 2 control points per sector): every leaf pair delivers its fluence profile with a unidirectional sliding window (`SlidingWindow`). Consecutive sectors sweep in opposite directions and meet at a shared leaf position, so the leaves move continuously along the arc. The slice time (the gantry slows down when needed), the sweep speed (machine leaf speed x slice time) and the MU per control point are solved together so that the sector needs the fewest MU: the MU of a leaf pair are its sum of upward fluence steps plus what is delivered while it travels, and the busiest row sets the MU of the sector. The sweep takes at most `maxSweepFraction` of the sector, and rows whose upward steps exceed `maxRiseRatio` x the sector peak are smoothed (mass preserving) first (`--rise-ratio`). The sampled trajectories are then refined (`fitIterations`) so the fluence delivered by the moving leaves matches the optimized fluence.
   - **LevelSet**: intensity-level decomposition, one segment per control point (fallback, and `--sequencer levels`).
   Control point apertures hold the leaf tips at the slice centres. With `integrateLeafMotion` (default, `--no-icm` to disable) the dose of a control point is the average over the leaf path through the neighbouring control points (`ApertureFluence::integratedOpenFractions`); without it the sweep cannot be represented, because a static aperture only reproduces about as many fluence levels as there are control points per sector.
3. `ArcPipeline::run` optimizes the aperture weights (MU) with L-BFGS (`weightIterations`) against a matrix-free or stored per-control-point dose operator, with MU smoothing between neighbouring control points and the MU penalty, then refines leaf positions and MU (step 5) and recomputes the dose directly.
4. `ArcSequencer::computeTiming` derives gantry speed, dose rate and leaf speed from the machine constraints and reports slowed control points.

5. Leaf/MU refinement (`--leaf-iter`, default 300): projected L-BFGS (memory 10, backtracking Armijo along the projection onto the speed/MU bounds, gradient restart on a bad direction) over all leaf tips and MU on the dose objectives, with exact tip gradients for the integrated motion. It replaced a Barzilai-Borwein gradient descent: on the lung case the sweep result went from PTV max 69.2 Gy / std 0.95 Gy / 770 MU to 68.3 Gy / 0.60 Gy / 374 MU. Tips move continuously inside the leaf-speed and MU-rate bounds of each slice and are snapped to the MLC resolution at the end.

### Progressive-resolution direct aperture optimization (`--strategy dao`)

This is the approach of Eclipse's VMAT optimizer (progressive resolution sampling, Otto 2008). It does not use the fluence sweep:

1. The fluence optimization still provides the target outline: a **conformal seed** (`ArcSequencer` method `Conformal`) opens every leaf pair over the extent of the sector fluence, at widely spaced control points (`daoCoarsestSpacingDeg`, default 16 deg).
2. **Leaf positions and MU are optimized directly** against the dose objectives (same exact-gradient matrix-free dose model, integrated leaf motion) with a bounded L-BFGS (`refineLeaves`).
3. Control points are **inserted by interpolation** (`ArcSequencer::resample`: tips linear in the arc angle, MU per degree linear) and the optimization continues. Spacings halve until the plan's control point spacing: 16, 8, 4, 2 deg. Coarse levels are cheap and allow large coordinated leaf moves; the finest level (twice the iterations) polishes.
4. The hot-spot term tightens from the fluence stage's setting to `daoHotspotThresholdPct` / `daoHotspotPenalty` (default 1.03 x Rx, 5000) over the resolutions. MU is not penalized by default (`daoMuPenalty`): open conformal apertures need little MU.

Leaf-speed and ordering limits are quadratic penalties during the iterations (box limits on tips and MU are clamped), and an exact projection repairs the result at the end of every resolution. A first version enforced the speed limit by projecting after each step; that turned descent directions into ascent directions and made most line searches fail (the objective then stalled at ~100x its final value), so it was replaced.

Lung case (HD120, 2 deg control points, 10 deg optimization beams, 500 iterations per resolution):

| | Fluence plan | Sweep | DAO |
|---|---|---|---|
| PTV D95 / D98 (Gy) | 66.1 / 66.0 | 65.8 / 65.5 | 66.0 / 65.3 |
| PTV max / std (Gy) | 68.7 / 0.37 | 68.4 / 0.66 | 70.7 / 1.14 |
| MU | | 471 | 95 |
| Delivery time | | 113 s, gantry slowed on 175 of 179 control points | 59.7 s, never slowed |
| Mean leaf travel per control point | | 18.8 mm (every leaf pair crosses the field each sector) | 1.4 mm |

The DAO gives the clinical look (open aperture around the target, small leaf motion, constant gantry speed, ~1/4 of the MU) but cannot reach the sweep's homogeneity: with one opening per leaf pair the modulation only comes from the angular dimension. A stricter hot-spot term lowers the maximum at the cost of coverage (1.01 x Rx / 50000: max 69.8 Gy, D95 63.7 Gy). Using only two resolutions (8, 4, 2) was worse (D95 63.7 Gy, std 2.1 Gy). A full DAO takes about 35 min (the finest resolution dominates; `--dao-iter 250` halves it at lower quality).

Hot spots: the delivered arc dose differs from the optimized fluence dose, and the difference shows up as hot spots. The arc stages therefore use their own, stricter hot-spot term (`hotspotThresholdPct` 1.01 x Rx, `hotspotPenalty` 50000; `--nto-threshold`, `--nto-penalty`) instead of the fluence stage's 1.04 / 2000. On the lung case this takes the PTV max from 71.6 to 69.2 Gy (fluence plan: 68.7) and the PTV std from 1.33 to 0.95 Gy with D95 65.5 Gy. Much stricter values (200000 above 1.00) start to cost coverage.

Conformal (experimental, `--sequencer conformal`): apertures follow the target outline and the leaf/MU optimization is left to add the modulation. The seed needs about a tenth of the MU of the sweep (87 vs 770) and no gantry slowdown, but the leaf optimization alone does not recover the fluence modulation (PTV std ~3 Gy, D95 ~60 Gy, also with L-BFGS, a 30x MU cap and 1000 iterations), so it is not the default. A real aperture-splitting DAO would be needed to make it competitive.

MU efficiency: a sliding window needs more MU than the fluence it delivers. `muPenalty` (`--mu-penalty`, default 0.3) adds a total-MU cost to the aperture weight and leaf stages, expressed as a fraction of the stage's initial objective so it does not depend on the objective scale. On the lung case it cuts the MU by about a third without worsening PTV D95. MU are in the units of the optimized fluence weights (course total), so the reported delivery time follows the machine's dose-rate limit applied to them.

Limits: DICOM export of VMAT plans is not implemented. The `leafSeq` CLI command refuses VMAT plans; use `arcSeq`. CLI: `plan --delivery vmat --arc <start> <stop> <cw|ccw> [couch] --opt-spacing <deg> --cp-spacing <deg>`.
