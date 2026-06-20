# `clinical.jl` - Clinical Genomics and Survival Analysis

## Overview

`clinical.jl` supports cohort containers, MAF parsing, TCGA file utilities, Kaplan-Meier and Cox survival analysis, ROC/CIF curves, pharmacogenomics helpers, and clinical visualization payloads.

### Purpose

This page documents the public surface of `clinical.jl`: the result types, workflow functions, and helper APIs a BioToolkit user is expected to call directly. Private helpers are intentionally omitted unless they define behavior users must understand.

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

## 1. Core Types

These containers hold clinical cohorts, mutation summaries, and model outputs.

| API | Description |
|---|---|
| `PatientCohort` | Patient/sample table with outcomes, covariates, and molecular annotations. |
| `KaplanMeierResult` | Survival curve estimates by group. |
| `CoxResult` | Cox model result with fitted terms and diagnostics. |
| `CoxTermResult` | One Cox coefficient, hazard ratio, interval, and p-value. |
| `MAFRecord` | One mutation annotation format record. |
| `MAFSummary` | Per-sample/per-gene mutation summary. |
| `ROCResult` | Sensitivity/specificity curve and AUC metadata. |
| `CIFResult` | Cumulative incidence function estimate. |
| `NeuralCoxResult` | Neural survival-model output. |
| `DoseResponseResult` | Dose-response fit and IC/EC-style summaries. |
| `OncoprintResult` | Matrix and annotation payload for oncoprints. |

## 2. Mutation and TCGA

MAF and TCGA helpers turn public cancer-genomics files into analysis-ready tables.

| API | Description |
|---|---|
| `read_maf` | Reads MAF records from a file or stream. |
| `summarize_maf` | Summarizes mutation burden, top genes, variant classes, and samples. |
| `tcga_query` | Builds a TCGA/GDC query description. |
| `tcga_download_files` | Downloads or stages files from a TCGA query result. |
| `merge_tcga_count_files` | Merges TCGA count files into a count matrix. |
| `tcga_ingest` | Runs query, download/staging, merge, and metadata assembly. |

## 3. Survival and Clinical Models

Model helpers return result objects and plotting payloads.

| API | Description |
|---|---|
| `kaplan_meier` | Computes Kaplan-Meier survival estimates. |
| `kaplan_meier_plot` | Builds a plot-ready survival curve payload. |
| `logrank_test` | Compares survival curves with a log-rank test. |
| `cox_ph` | Fits Cox proportional hazards models. |
| `forest_plot` | Creates forest-plot data from Cox or other interval results. |
| `survival_roc` | Computes time-dependent survival ROC summaries. |
| `cif_curve` | Estimates cumulative incidence under competing risks. |
| `neural_cox` | Fits a compact neural Cox-style model. |
| `dose_response_curve` | Fits and summarizes dose-response data. |
| `oncoprint` | Builds an alteration matrix for cohort-level mutation display. |

## 4. Pharmacogenomics and Cohorts

These helpers support translational cohort interpretation.

| API | Description |
|---|---|
| `pharmacogenomics_star_alleles` | Calls or summarizes star-allele diplotypes from variant inputs. |
| `cpic_metabolizer_phenotype` | Maps diplotypes to CPIC-like metabolizer phenotypes. |
| `pharmgkb_like_recommendations` | Creates PharmGKB-style drug/gene recommendation summaries. |
| `trial_suitability_scores` | Scores patient records against trial inclusion/exclusion features. |
| `propensity_score_match` | Performs simple propensity-score matching. |
| `omop_visit_summary` | Summarizes OMOP visit records. |
| `synthpop_like_cohort` | Generates synthetic cohort-like tables for testing workflows. |

---

## Complete Usage Example

```julia
using BioToolkit

maf = read_maf("cohort.maf")
summary = summarize_maf(maf)
km = kaplan_meier(times, events, groups)
fit = cox_ph(cohort, :survival_time, :event; covariates=[:age, :stage])
```

