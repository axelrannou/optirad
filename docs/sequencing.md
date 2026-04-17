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
