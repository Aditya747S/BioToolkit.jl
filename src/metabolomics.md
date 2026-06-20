# `metabolomics.jl` - Metabolomics Analysis

## Overview

`metabolomics.jl` implements metabolomics normalization, adduct and isotope annotation, pathway enrichment, QC, PCA, correlation networks, source tracking, differential abundance, NMR peak tables, targeted quantification, and power analysis.

### Purpose

This page documents the public API in `metabolomics.jl`, grouped by the biological workflow it supports. Each entry describes the role of the function or type in practical BioToolkit analyses.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Domain-specific names** | Public functions match familiar analysis concepts. |
| **Compact result containers** | Typed results preserve the fields needed for downstream inspection. |
| **Composable tables and matrices** | APIs accept common Julia data structures and BioToolkit containers. |
| **Deterministic summaries** | Statistical summaries are reproducible unless stochastic behavior is explicitly requested. |
| **Workflow interoperability** | Outputs can feed plotting, enrichment, clinical, or systems biology modules. |

---

## 1. Types

Typed results preserve source tracking, feature annotation, isotope tracing, and pathway enrichment outputs.

| API | Description |
|---|---|
| `MetabolomicsSourceTrackingResult` | Posterior/source-contribution summary for metabolomics source tracking. |
| `MetaboliteAnnotation` | Candidate metabolite identity with mass/adduct/error metadata. |
| `IsotopeTrace` | Isotope tracing result for one feature or pathway. |
| `MetabolicPathwayResult` | Pathway enrichment/activity result from metabolite sets. |

## 2. Preprocessing and Annotation

These functions prepare feature matrices and annotate metabolite identities.

| API | Description |
|---|---|
| `normalize_metabolomics` | Normalizes intensity matrices by TIC, median, quantile, or supplied factors. |
| `quality_control_metabolomics` | Computes QC metrics and flags suspect samples/features. |
| `missing_value_analysis` | Summarizes missingness and missingness mechanisms. |
| `batch_effect_correction_metabolomics` | Corrects batch effects in metabolomics matrices. |
| `metabolomics_streaming_analysis` | Processes feature records in streaming mode. |
| `annotate_metabolite_features` | Annotates features against mass/name databases. |
| `annotate_adducts` | Annotates likely adduct forms. |
| `isotope_tracer_analysis` | Analyzes isotope-label incorporation. |
| `metabolite_classification` | Assigns metabolite classes from annotation fields. |
| `nmr_peak_table` | Builds NMR peak table summaries. |
| `targeted_quantification` | Quantifies targeted metabolites from feature tables. |

## 3. Statistics and Pathways

Analysis functions support differential abundance, pathway interpretation, and network summaries.

| API | Description |
|---|---|
| `metabolomics_source_tracking_model` | Fits a source-tracking model. |
| `metabolomics_source_tracking` | Estimates source contributions. |
| `metabolomics_source_tracking_posterior_summary` | Summarizes source-tracking posterior output. |
| `metabolomics_differential_abundance` | Tests metabolite differential abundance between groups. |
| `quantify_metabolite_variation` | Summarizes variance components or feature variability. |
| `metabolite_pathway_enrichment` | Runs pathway enrichment over annotated metabolites. |
| `pca_metabolomics` | Computes PCA for metabolomics matrices. |
| `metabolite_correlation_network` | Builds metabolite correlation networks. |
| `fold_change_metabolomics` | Computes group fold changes. |
| `volcano_metabolomics` | Creates volcano-plot-ready metabolomics results. |
| `metabolomics_power_analysis` | Estimates power/sample size for metabolomics contrasts. |
| `feature_clustering_metabolomics` | Clusters metabolite features. |

---

## Complete Usage Example

```julia
using BioToolkit

norm = normalize_metabolomics(intensity_matrix)
qc = quality_control_metabolomics(norm)
ann = annotate_adducts(feature_table)
da = metabolomics_differential_abundance(norm, group_labels)
```

