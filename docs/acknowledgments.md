# Acknowledgments and Upstream References

OptiRad was developed as an independent C++ radiotherapy treatment planning system, but several algorithms, data structures, and implementation conventions were derived from or informed by matRad.

## matRad

matRad is an open-source radiotherapy treatment planning system developed by the Radiotherapy Optimization Group at the German Cancer Research Center (DKFZ).

- Project: matRad
- Upstream: https://github.com/e0404/matRad
- License: BSD 3-Clause License
- Copyright: the matRad development team

Portions of OptiRad were implemented by translating, adapting, or closely following concepts from matRad's MATLAB code and data model. This document records those relationships to make the provenance of major planning components explicit.

## Directly Ported or Closely Adapted Components

### Dose Calculation

| OptiRad component | matRad reference | Notes |
|---|---|---|
| `SSDCalculator` | `matRad_computeSSD.m` | Source-to-surface distance calculation logic, including nearest valid-ray fallback behavior. |
| `RadDepthCalculator` | `matRad_rayTracing.m` | Radiological depth calculation along beam paths. |
| `PencilBeamEngine` | `matRad_PhotonPencilBeamSVDEngine` | Photon pencil-beam dose engine using matRad-aligned kernel conventions, lateral cutoffs, and interpolation layout. |
| `DoseCalcOptions` | matRad batched dose workflows | Thresholding and batch-style calculation options were inspired by matRad's dose pipeline. |

### Plan Analysis and Optimization

| OptiRad component | matRad reference | Notes |
|---|---|---|
| `PlanAnalysis` | `matRad_planAnalysis.m` | DVH curves, Dx%, VxGy, conformity index, and homogeneity index reporting. |
| `DVHObjective` | matRad DVH penalty formulation | Voxel-sorting DVH penalty structure follows the same formulation approach. |

### Sequencing and Deliverability

| OptiRad component | matRad reference | Notes |
|---|---|---|
| `DeliverableDoseCalculator` | matRad sequencing workflow | Deliverable weight reconstruction follows the same quantized-fluence interpretation used by matRad. |
| Step-and-shoot sequencing flow | `matRad_sequencing.m`, `matRad_directApertureOptimization.m`, `matRad_recalcApertureInfo.m` | Overall workflow and terminology are aligned with matRad sequencing concepts. |

## Data Structures and Conventions Inspired by matRad

| OptiRad component | matRad equivalent | Notes |
|---|---|---|
| `PatientData` | `ct` + `cst` | Central patient container mirrors the separation of image data and structure definitions. |
| `Structure` voxel indices | `cst{i,4}` | Precomputed voxel membership follows matRad's structure representation pattern. |
| `Structure` optimization parameters | `cst{i,5}` | Per-structure optimization settings follow the same conceptual split. |
| `StfProperties` | `pln.propStf` | Beam-angle and steering configuration follows matRad-style treatment steering properties. |
| `PhotonIMRTStfGenerator` | `matRad_generateStf.m` | Beam generation behavior and explicit gantry/couch list support were influenced by matRad. |
| Grid interpolation behavior | matRad interpolation conventions | Several nearest-neighbor and grid-indexing behaviors intentionally match matRad semantics. |
| Rotation and coordinate conventions | matRad geometry utilities | Some transformation conventions were chosen to ease comparison with matRad results. |

## Why This Attribution Exists

This project contains comments and design decisions that explicitly reference matRad. Recording that relationship is useful for:

- scientific transparency,
- reproducibility,
- license compliance,
- and helping future contributors understand why certain behaviors intentionally match matRad.

## License Notice for Derived Material

matRad is distributed under the BSD 3-Clause License. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the upstream copyright notice, license conditions, and disclaimer are preserved.

When redistributing OptiRad source or binaries, retain the matRad attribution described in this document and preserve upstream notices where required.

For the authoritative upstream license text, see:
https://github.com/e0404/matRad/blob/master/LICENSE.md