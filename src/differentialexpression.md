# `differentialexpression.jl` - Differential Expression Analysis

## Overview

`differentialexpression.jl` implements DESeq2-, edgeR-, and limma-inspired workflows for count matrices: normalization, dispersion estimation, negative-binomial testing, GLM/QL tests, effect-size shrinkage, transformations, and batch correction.

### Purpose

This file documents the exported analysis objects and workflow functions in `differentialexpression.jl`. It focuses on what each API does, what stage of the workflow it belongs to, and how the pieces compose with the rest of BioToolkit.

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

## 1. Core Types and Setup

Dataset containers preserve counts, design matrices, feature names, sample names, and model state.

| API | Description |
|---|---|
| `CountMatrix` | Count matrix with gene/sample identifiers and metadata. |
| `DEResult` | Differential-expression row/result object. |
| `GLMSolver` | Configuration for GLM fitting routines. |
| `DESeqDataSet` | DESeq-style dataset with counts, design, size factors, dispersions, and results. |
| `DESeqDataSetFromMatrix` | Builds a `DESeqDataSet` from count matrix and sample metadata/design. |
| `counts` | Retrieves count data from a dataset. |
| `design` | Retrieves design information. |
| `makeExampleDESeqDataSet` | Creates a small synthetic DESeq-style dataset for examples/tests. |

## 2. Normalization and Dispersion

These functions estimate size factors, normalization factors, and dispersion trends/priors.

| API | Description |
|---|---|
| `calc_tmm_factors` | Computes TMM-like library normalization factors. |
| `calc_norm_factors` | General normalization-factor wrapper. |
| `estimateSizeFactorsForMatrix` | Median-ratio size factor estimation from a raw matrix. |
| `estimateSizeFactors` | Stores size factors in a dataset. |
| `sizeFactors` | Reads size factors. |
| `sizeFactors!` | Sets size factors. |
| `normalizationFactors` | Reads gene/sample normalization factors. |
| `normalizationFactors!` | Sets normalization factors. |
| `estimate_dispersions` | Estimates per-gene negative-binomial dispersions. |
| `estimate_dispersions_prior` | Estimates prior/trend information for dispersion shrinkage. |
| `estimateDispersions` | Full DESeq-style dispersion workflow. |
| `estimateDispersionsGeneEst` | Gene-wise dispersion estimation. |
| `estimateDispersionsFit` | Fits a dispersion trend. |
| `estimateDispersionsMAP` | MAP/shrunken dispersion estimation. |
| `estimateDispersionsPriorVar` | Estimates dispersion prior variance. |

## 3. Testing and Results

Wald, likelihood-ratio, exact, and quasi-likelihood workflows are available.

| API | Description |
|---|---|
| `nbinomWaldTest` | Runs negative-binomial Wald tests. |
| `nbinomLRT` | Runs negative-binomial likelihood-ratio tests. |
| `DESeq` | Runs the standard DESeq-style pipeline. |
| `results` | Extracts contrast/coefficient results. |
| `resultsNames` | Lists available result coefficients. |
| `benjamini_hochberg` | Adjusts p-values by BH FDR. |
| `differential_expression` | Convenience wrapper for count filtering, normalization, testing, and results. |
| `filter_low_counts` | Removes weakly expressed genes. |
| `shrink_lfc` | Applies log-fold-change shrinkage. |
| `cooks_distance` | Computes Cook distance diagnostics. |
| `replace_outliers` | Replaces flagged count outliers. |

## 4. edgeR, Transforms, and Batch Effects

Additional workflows support edgeR-like tests and downstream matrix correction.

| API | Description |
|---|---|
| `estimate_dispersion_edgeR` | Estimates edgeR-style dispersions. |
| `trend_dispersion_edgeR` | Fits dispersion trends. |
| `dispersion_prior_dof` | Estimates prior degrees of freedom. |
| `exact_test_edgeR` | Runs exact tests for two-group count comparisons. |
| `glm_ql_fit` | Fits quasi-likelihood GLMs. |
| `glm_ql_f_test` | Runs QL F-tests. |
| `edgeR_qlf_test` | Convenience edgeR QL workflow. |
| `vst` | Variance-stabilizing transformation. |
| `fit_gene_fast!` | Fits one gene/model efficiently in place. |
| `remove_batch_effect` | Linear batch-effect removal. |
| `combat_correction` | ComBat-style empirical Bayes batch correction. |
| `estimate_surrogates` | Estimates surrogate variables from expression matrices. |

---

## Complete Usage Example

```julia
using BioToolkit

dds = DESeqDataSetFromMatrix(counts_matrix, sample_table, design_matrix)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
res = results(dds)
```

