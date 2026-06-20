# `epigenetics.jl` - Epigenomics Analysis

## Overview

`epigenetics.jl` covers coverage vectors, peak calling, differential binding, methylation analysis, single-cell chromatin transforms, coaccessibility, motif deviations, footprints, and Hi-C compartment/TAD helpers.

### Purpose

This file documents the exported analysis objects and workflow functions in `epigenetics.jl`. It focuses on what each API does, what stage of the workflow it belongs to, and how the pieces compose with the rest of BioToolkit.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **End-to-end workflow coverage** | Public functions cover preprocessing, modeling, diagnostics, and reporting. |
| **Explicit result objects** | Important outputs are typed so fields can be inspected and reused. |
| **Method-compatible inputs** | Matrix, table, and BioToolkit container inputs are supported where the operation naturally allows it. |
| **Statistical transparency** | Helpers expose normalization, filtering, testing, and correction steps rather than hiding them. |
| **Plot/report readiness** | Returned objects are structured for downstream plotting and export. |

---

## 1. Core Types

Typed containers represent genomic signal, peaks, methylation, chromatin experiments, and contact matrices.

| API | Description |
|---|---|
| `SparseCoverageVector` | Sparse genomic coverage representation. |
| `Peak` | Genomic peak interval with score/statistics. |
| `PeakSet` | Collection of peaks with metadata. |
| `PeakSupport` | Support evidence for a peak across samples/replicates. |
| `Epigenome` | Combined epigenomic tracks and annotations. |
| `MethylationCall` | One methylation observation. |
| `MethylationExperiment` | Methylation matrix plus sample/feature metadata. |
| `MethylationResult` | Differential methylation result. |
| `SingleCellChromatinExperiment` | Single-cell accessibility/chromatin container. |
| `CoaccessibilityEdge` | Coaccessible peak pair with score. |
| `TadResult` | Detected TAD interval and score. |
| `ContactMatrix` | Hi-C/contact matrix with bin metadata. |

## 2. Coverage, Peaks, and Binding

Signal-processing functions convert alignments/intervals into coverage and peak-level tests.

| API | Description |
|---|---|
| `calculate_coverage` | Computes coverage from intervals or aligned records. |
| `coverage_depth` | Summarizes depth over regions. |
| `coverage_segments` | Converts coverage to contiguous segments. |
| `call_peaks` | Calls enriched intervals from coverage/signal. |
| `normalize_gc_bias` | Corrects signal for GC-content bias. |
| `summarize_peak_support` | Summarizes peak reproducibility/support. |
| `count_overlaps` | Counts overlaps between peaks and reads/features. |
| `differential_binding` | Tests peak-level differential accessibility/binding. |

## 3. Methylation and Single-Cell Chromatin

Methylation and scATAC helpers support common preprocessing and tests.

| API | Description |
|---|---|
| `bin_methylation` | Aggregates methylation calls into genomic bins. |
| `differential_methylation` | Tests methylation differences between groups. |
| `parse_bismark_coverage` | Reads Bismark coverage-style methylation files. |
| `tfidf` | Applies TF-IDF normalization to accessibility counts. |
| `run_lsi` | Runs latent semantic indexing. |
| `rsvd` | Randomized SVD helper used by LSI. |
| `gene_activity_score` | Computes gene activity from nearby accessibility. |
| `calculate_coaccessibility` | Scores peak coaccessibility. |
| `compute_motif_deviations` | Computes motif accessibility deviation scores. |
| `detect_footprints` | Detects motif footprint patterns. |
| `chromhmm_like_segmentation` | Segments chromatin states from signal tracks. |
| `tobias_like_footprints` | Computes TOBIAS-style footprint summaries. |

## 4. Hi-C and 3D Genome

Contact-matrix helpers normalize maps and call domains/compartments.

| API | Description |
|---|---|
| `directionality_index` | Computes directionality index across contact bins. |
| `detect_tads` | Calls topologically associating domains. |
| `hic_ice_normalize` | Performs ICE-style matrix balancing. |
| `hic_kr_normalize` | Performs KR-style matrix balancing. |
| `hic_ab_compartments` | Computes A/B compartment scores. |
| `insulation_score` | Computes insulation scores across windows. |
| `compartment_boundary_candidates` | Finds candidate compartment boundaries. |

---

## Complete Usage Example

```julia
using BioToolkit

coverage = calculate_coverage(intervals)
peaks = call_peaks(coverage)
counts = count_overlaps(peaks, fragments)
lsi = run_lsi(tfidf(counts))
```

