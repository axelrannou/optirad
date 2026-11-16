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

## VMAT: direct aperture optimization (DAO)

Delivery mode is chosen at planning time (`SMLC` or `VMAT`, `Plan::setDeliveryMode`). For VMAT the arc is defined by start angle, stop angle, direction and couch angle (`ArcSpec`), plus an arc sector width (default 20 deg; one STF beam per sector, used for geometry only) and a control-point spacing (default 4 deg) (`VmatSettings`).

VMAT plans are always optimized with progressive-resolution DAO, the approach of Eclipse's VMAT optimizer (progressive resolution sampling, Otto 2008). No Dij and no fluence optimization are needed: DAO runs straight after planning.

1. **Conformal seed** (`ArcSequencer::conformalSeed`): every leaf pair is opened over the target projection at widely spaced control points (`coarsestSpacingDeg`, default 16 deg) with 1 MU each; the MU scale is then calibrated against the objectives.
2. **Leaf positions and MU are optimized directly** against the dose objectives and the NTO with a bounded L-BFGS (`refineLeaves`), using an exact-gradient matrix-free dose model with integrated leaf motion (`integrateLeafMotion`, `--no-icm` to disable).
3. **Control points are inserted by interpolation** (`ArcSequencer::resample`: tips linear in the arc angle, MU per degree linear) and the optimization continues. Spacings halve until the plan's control point spacing (e.g. 16, 8, 4, 2 deg). Coarse levels allow large coordinated leaf moves; the finest level (twice the iterations) polishes.
4. `ArcSequencer::computeTiming` derives gantry speed, dose rate and leaf speed from the machine constraints and reports slowed control points.

**MU units.** Control point `mu` is a bixel weight for the whole course (the objectives are in total course dose). The machine's dose-rate limits (MU/s per fraction) are converted with the matRad calibration of the photon base data, 1 weight unit per fraction = 100 MU × (SAD / 1000 mm)², i.e. `Plan::machineMuPerWeight()` (also `ArcPipeline::machineMuPerWeight(plan)`) = 100 (SAD/1 m)² / fractions (3.33 for 30 fractions at SAD 1000 mm), and `ArcPipeline::deliveryConstraints(plan)` holds the converted limits used by the DAO bounds and the timing. The pipeline logs the total as approximate machine MU per fraction. The upper MU bound of a control point lets the gantry slow down (to the machine minimum gantry speed, or 1/10 of the maximum when the minimum is 0); leaf travel limits keep the full-speed time.

The L-BFGS direction is computed on the projected gradient (variables on an active bound are excluded), and each resolution stops at its iteration budget or when the best cost improves less than 0.1 % over 40 accepted steps.

Leaf-speed and ordering limits are quadratic penalties during the iterations (box limits on tips and MU are clamped), and an exact projection repairs the result at the end of every resolution (projecting after each step turned descent directions into ascent directions). MU smoothing between neighbouring control points (`muSmoothingWeight`) and an optional total-MU penalty (`muPenalty`, fraction of the initial objective, default 0) are added to the cost.

The selected compute backend also applies to the final delivered-dose calculation. `Auto` tries CUDA and logs a warning if it falls back to CPU; `CUDA` requires GPU execution; `CPU` disables it.

Hot spots: there is no global hot-spot penalty. As in Eclipse, control target hot spots with Upper / Max DVH objectives on the PTV and the dose fall-off outside the target with the normal tissue objective (see `docs/optimization.md`).

With one opening per leaf pair, modulation comes from the angular dimension only: DAO gives the clinical look (open aperture around the target, small leaf motion, constant gantry speed, low MU) at the cost of some target homogeneity compared with an ideal fluence plan.

Limits: DICOM RT Plan export of VMAT plans is not implemented (the IMRT step-and-shoot export writes `BeamMeterset` in machine MU per fraction with the same calibration; control point cumulative weights stay relative). The `leafSeq` CLI command refuses VMAT plans; use `arcSeq`. CLI: `plan --delivery vmat --arc <start> <stop> <cw|ccw> [couch] --opt-spacing <sector deg> --cp-spacing <deg>`, then `arcSeq [--dose-resolution <mm>] [--backend ...] [--dao-iter <n>] [--dao-coarsest <deg>] [--mu-penalty <f>] [--mu-smooth <w>] [--nto-* ...]`.
