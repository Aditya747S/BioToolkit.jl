# systemsbio.jl

## Purpose
This file implements systems biology routines for network inference, module detection, limma-style modeling, gene-set enrichment, and multi-omics factor analysis.

## Main types
- `GeneNetwork` stores a weighted gene graph, node mappings, connectivity scores, and module assignments.
- `SoftThresholdResult` stores the power sweep used for scale-free network selection.
- `LimmaResult` stores coefficients, errors, statistics, p-values, q-values, variances, and labels.
- `NetworkInferenceResult` stores a directed graph and BIC score for inferred networks.
- `MultiOmicsFactorAnalysisResult` stores factor loadings and explained variance across assays.

## Public functions
- `pick_soft_threshold(counts; powers)` evaluates scale-free fit across candidate powers.
- `find_modules(counts; power, tom_threshold, edge_threshold, min_module_size)` builds a gene co-expression network and module assignment.
- `module_eigengenes(counts, modules)` computes module-level eigengenes.
- `module_dendrogram(network)` converts module assignments into a tree-like summary.
- `limma_fit(counts, design; coefficient_names)` fits a moderated linear model.
- `limma_deresults(result; coefficient_index)` converts a `LimmaResult` into `DEResult` rows.
- `gsea(results, database; permutations, min_size, max_size, seed)` runs gene-set enrichment on ranked results.
- `infer_network` and `multi_omics_factor_analysis` are exported analysis entry points for network inference and multi-omics integration.

## How it is used
The usual workflow is to start with count data, call `pick_soft_threshold` or `find_modules` to build a co-expression network, and then use `module_eigengenes` or `module_dendrogram` to summarize the module structure.

For differential modeling, `limma_fit` and `limma_deresults` provide moderated linear statistics that can be passed straight into `gsea`. The network and factor-analysis entry points let the same module cover both graph-based and latent-factor systems biology workflows.

## Implementation notes
- The network routines depend on weighted and unweighted graph primitives from `Graphs` and `SimpleWeightedGraphs`.
- `limma_fit` uses empirical Bayes-style variance moderation.
- `gsea` applies Benjamini-Hochberg correction after permutation-based p-value estimation.

## Why it matters
Systems biology is where single-gene measurements become interpretable modules and pathways. This file provides the network, model, and enrichment pieces needed to make that transition inside BioToolkit.
