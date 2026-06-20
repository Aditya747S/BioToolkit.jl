# `somatic.jl` - Somatic Variant Analysis

## Overview

`somatic.jl` provides germline and tumor-normal calling helpers, external caller wrappers, structural-variant breakpoint summaries, VEP-like annotation, tumor mutational burden, mutational signatures, copy-number segmentation, clonal deconvolution, driver enrichment, variant tiering, hotspot scans, and mutational spectra.

### Purpose

This documentation covers the public API implemented in `somatic.jl`, with concrete descriptions for the result types and workflow functions exposed by BioToolkit.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Concrete workflow coverage** | Each section maps to a real analysis task implemented by the module. |
| **Typed results** | Important model outputs and summaries are represented with explicit structs. |
| **Downstream compatibility** | Results are shaped for plotting, tabulation, or reuse in other BioToolkit modules. |
| **Source-aligned API names** | Entries use implemented/exported names rather than speculative aliases. |
| **Readable defaults** | Examples show the minimal flow without hiding required biological inputs. |

---

## 1. Calling and Wrappers

These functions produce or prepare variant calls from germline, tumor-normal, and external caller workflows.

| API | Description |
|---|---|
| `germline_bayesian_call` | Bayesian germline genotype/variant calling summary. |
| `somatic_tumor_normal_call` | Tumor-normal somatic variant calling summary. |
| `haplotypecaller_call` | Wrapper/payload helper for HaplotypeCaller-style calls. |
| `mutect2_call` | Wrapper/payload helper for Mutect2-style calls. |
| `strelka2_call` | Wrapper/payload helper for Strelka2-style calls. |

## 2. Annotation, Burden, and Tiering

Variant interpretation helpers annotate and prioritize somatic events.

| API | Description |
|---|---|
| `vep_like_annotation` | Adds consequence/gene/impact annotations in a VEP-like format. |
| `tumor_mutational_burden` | Computes mutations per callable megabase. |
| `variant_tier_classification` | Assigns clinical/research tiers from consequence and evidence fields. |
| `somatic_hotspot_scan` | Finds recurrent hotspot mutations. |
| `driver_enrichment` | Tests driver-gene or pathway enrichment. |
| `mutational_spectrum_96` | Builds 96-channel trinucleotide mutational spectrum. |

## 3. Signatures, Copy Number, and Clonality

Higher-level tumor evolution functions summarize mutational processes and clones.

| API | Description |
|---|---|
| `mutational_signature_nmf` | Fits NMF signatures to mutation spectra. |
| `cosmic_signature_attribution` | Attributes spectra to COSMIC-like signatures. |
| `cbs_like_segments` | Segments copy-number signal with CBS-like behavior. |
| `copy_number_from_coverage` | Estimates copy number from coverage depth. |
| `allelic_imbalance_test` | Tests allele imbalance from counts or BAF values. |
| `sv_breakpoint_graph_table` | Builds structural-variant breakpoint graph table. |
| `sv_fusion_candidates` | Finds candidate gene fusions from SV breakpoints. |
| `clonal_deconvolution_pyclone` | PyClone-like clonal deconvolution summary. |
| `clonal_evolution_tree` | Builds a clonal evolution tree. |

---

## Complete Usage Example

```julia
using BioToolkit

calls = somatic_tumor_normal_call(tumor_counts, normal_counts)
ann = vep_like_annotation(calls, gene_features)
tmb = tumor_mutational_burden(ann; callable_mb=35.0)
sig = mutational_signature_nmf(mutational_spectrum_96(ann))
```

