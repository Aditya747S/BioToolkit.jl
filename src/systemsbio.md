# `systemsbio.jl` - Systems Biology and Multi-Omics

## Overview

`systemsbio.jl` implements WGCNA-like module discovery, topological overlap, module eigengenes, limma/voom modeling, empirical Bayes moderation, duplicate correlation, batch-effect removal, GSEA integration, Bayesian network inference, MOFA/CCA-style multi-omics integration, WNN-like integration, GRN inference, NMF gene programs, sparse CCA, and causal regression.

### Purpose

This documentation covers the public API implemented in `systemsbio.jl`, with concrete descriptions for the result types and workflow functions exposed by BioToolkit.

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

## 1. Core Types

Result containers capture networks, limma/voom fits, inferred graphs, and multi-omics factors.

| API | Description |
|---|---|
| `GeneNetwork` | Weighted gene network with graph, gene mappings, modules, and connectivity. |
| `SoftThresholdResult` | Power scan result for scale-free topology and connectivity. |
| `LimmaResult` | Moderated linear-model coefficients, statistics, p-values, q-values, and metadata. |
| `VoomResult` | Log-CPM matrix, precision weights, fitted trend, effective library sizes, and design. |
| `NetworkInferenceResult` | Inferred directed network plus gene/node mappings and BIC score. |
| `MultiOmicsFactorAnalysisResult` | Latent factors, feature loadings, assay-specific loadings, and explained variance. |

## 2. Networks and Modules

Network functions build coexpression modules and graph summaries.

| API | Description |
|---|---|
| `pick_soft_threshold` | Scans powers to choose scale-free coexpression threshold. |
| `find_modules` | Builds adjacency/TOM network and assigns modules. |
| `module_eigengenes` | Computes first-PC eigengene per module. |
| `module_dendrogram` | Creates a phylogenetic-style module dendrogram. |
| `infer_network` | Infers directed regulatory/dependency network by score search. |
| `causal_grn_inference` | Builds causal-GRN candidate edges from expression associations. |
| `gene_program_nmf` | Finds nonnegative gene programs. |

## 3. limma, voom, and Enrichment

Statistical modeling functions implement limma-style empirical Bayes workflows.

| API | Description |
|---|---|
| `voom_transform` | Transforms counts to log-CPM and precision weights. |
| `voom_with_quality_weights` | Adds sample quality weights to voom output. |
| `estimate_limma_hyperparameters` | Estimates prior variance/degrees of freedom for moderation. |
| `limma_moderate_ttest` | Runs moderated t-test for one response/design. |
| `limma_fit` | Fits limma/voom model over a count matrix. |
| `limma_deresults` | Converts `LimmaResult` to differential-expression rows. |
| `eBayes` | Applies empirical Bayes moderation to a fit. |
| `contrasts_fit` | Applies contrast matrix to limma coefficients. |
| `duplicateCorrelation` | Estimates within-block duplicate correlation. |
| `remove_batch_effect_limma` | Removes batch effects with limma-style linear modeling. |
| `gsea` | Runs GSEA over DE or limma results. |

## 4. Multi-Omics Integration and Causal Models

Integration helpers compute shared latent factors, anchors, and causal estimates.

| API | Description |
|---|---|
| `multi_omics_factor_analysis` | Computes shared factors across assays. |
| `mofa_plus_integration` | MOFA+-style integration wrapper. |
| `mofa_plus_em` | Expectation-maximization MOFA-like factorization. |
| `cca_integration` | Canonical correlation analysis integration. |
| `sparse_cca_integration` | Sparse CCA integration. |
| `anchor_based_integration` | Finds anchors between reference and query assays. |
| `totalvi_like_integration` | Integrates RNA and protein counts in totalVI-like latent space. |
| `causal_ate_regression` | Estimates average treatment effect by regression adjustment. |

---

## Complete Usage Example

```julia
using BioToolkit

soft = pick_soft_threshold(counts)
network = find_modules(counts; power=soft.best_power)
voom = voom_transform(counts, design_matrix)
fit = limma_fit(counts, design_matrix)
factors = multi_omics_factor_analysis([rna_matrix, protein_matrix])
```

