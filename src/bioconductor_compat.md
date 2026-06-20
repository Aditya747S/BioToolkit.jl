# `bioconductor_compat.jl` - Bioconductor Compatibility Helpers

## Overview

`bioconductor_compat.jl` provides Julia-native equivalents or payload builders for common Bioconductor workflows across single-cell analysis, epigenomics, metabolomics, variants, GWAS, and visualization.

### Purpose

This page documents the public surface of `bioconductor_compat.jl`: the result types, workflow functions, and helper APIs a BioToolkit user is expected to call directly. Private helpers are intentionally omitted unless they define behavior users must understand.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Structured domain outputs** | Stable results are returned as structs, named tuples, or table-friendly payloads. |
| **Composable Julia inputs** | APIs favor ordinary arrays, matrices, dictionaries, BioToolkit records, and DataFrames. |
| **Workflow granularity** | Each public function maps to a recognizable analysis or data-preparation step. |
| **Optional external integration** | Wrapper-style functions prepare command inputs or parse results without making all workflows depend on external tools. |
| **Reporting-ready summaries** | Many functions return data that can be plotted, exported, or passed into downstream modules. |

---

## 1. Parallel and Single-Cell

Helpers mirror familiar Bioconductor package roles while returning Julia tables or matrices.

| API | Description |
|---|---|
| `biocparallel_map` | Applies a function over inputs with a BiocParallel-like interface. |
| `scater_qc_metrics` | Computes per-cell and per-feature quality-control summaries. |
| `scran_pool_size_factors` | Estimates pooled size factors in the style of scran. |
| `scuttle_aggregate_counts` | Aggregates counts by groups for pseudobulk or summary analysis. |
| `singler_annotate` | Assigns labels by reference expression similarity. |
| `scmap_project` | Projects query cells onto a reference embedding or feature space. |
| `splatter_simulate_counts` | Simulates count matrices with Splatter-like mean/dispersion structure. |

## 2. Epigenomics

These functions approximate common DiffBind, csaw, chromVAR, minfi, and bsseq workflows.

| API | Description |
|---|---|
| `diffbind_differential_binding` | Builds differential binding summaries from peak-count matrices and design labels. |
| `csaw_window_differential_binding` | Tests fixed genomic windows for differential signal. |
| `chromvar_deviations` | Computes motif deviation scores from accessibility and motif incidence matrices. |
| `minfi_dmp` | Finds differentially methylated positions from methylation beta/M-value matrices. |
| `bsseq_dmr` | Identifies region-level methylation differences. |

## 3. Metabolomics, Variants, and Plots

Compatibility outputs are shaped for downstream BioToolkit plotting and reporting.

| API | Description |
|---|---|
| `xcms_peak_workflow` | Groups and summarizes LC-MS peak tables in an xcms-like workflow. |
| `lipidr_differential_abundance` | Runs lipid/metabolite differential abundance summaries. |
| `mixomics_factor_analysis` | Computes mixOmics-style latent factors across assays. |
| `mofa2_factor_analysis` | Builds MOFA-like low-dimensional factors from multi-omics matrices. |
| `variantannotation_annotate` | Annotates variant records with feature overlaps and predicted consequences. |
| `varianttools_filter_variants` | Filters variants by quality, consequence, frequency, or supplied predicates. |
| `gwascat_lookup` | Creates GWAS Catalog-like lookup summaries. |
| `complexheatmap_payload` | Returns row/column ordering and annotation data for heatmap rendering. |
| `gviz_track_table` | Builds genomic track tables compatible with Gviz-style visualization. |

---

## Complete Usage Example

```julia
using BioToolkit

qc = scater_qc_metrics(counts_matrix)
size_factors = scran_pool_size_factors(counts_matrix)
labels = singler_annotate(query_expression, reference_expression, reference_labels)
```

