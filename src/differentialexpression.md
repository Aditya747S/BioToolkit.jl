# differentialexpression.jl

## Purpose
This file implements the count-based differential expression workflow in BioToolkit. It provides the core count-matrix container, normalization factor estimation, dispersion modeling, GLM fitting helpers, and result types used for RNA-seq style analysis.

## Main structs
- CountMatrix: a sparse gene-by-sample count matrix with gene IDs and sample IDs.
- DEResult: one differential expression result containing base mean, log2 fold change, standard error, test statistic, p-value, and adjusted p-value.
- GLMSolver: reusable internal state for iterative GLM fitting.

## Public functions
- CountMatrix(counts, gene_ids, sample_ids): construct a validated count matrix.
- filter_low_counts(counts; min_total): remove genes with very low total counts.
- calc_norm_factors(counts): estimate sample normalization factors.
- estimate_dispersions(counts, norm_factors): estimate gene-wise dispersion values.
- differential_expression: perform the full differential expression analysis.
- benjamini_hochberg: adjust p-values for multiple testing.
- shrink_lfc: shrink fold-change estimates.
- vst: compute a variance-stabilizing transformation.
- fit_gene_fast!: fit a gene-level GLM in place using a fast iterative solver.
- remove_batch_effect and combat_correction: correct batch effects.
- estimate_surrogates: estimate surrogate variables.

## What the module does
The workflow starts with raw count data and turns it into normalized, model-ready statistics. Low-count genes can be filtered first, then normalization factors and dispersions are estimated, and finally gene-wise GLM models are fitted to obtain test statistics and p-values. The output is a collection of DEResult records that can be converted to tables or fed into plotting and enrichment modules.

## How the structs work together
CountMatrix is the main input container and is used throughout the pipeline. GLMSolver is the iterative fitting state that stores design matrices, weights, fitted coefficients, and working variables. DEResult is the final user-facing output structure that captures the statistical evidence for each gene.

## Typical usage
1. Build a CountMatrix from a sparse integer matrix or a dense matrix.
2. Run filter_low_counts to remove weakly expressed genes.
3. Estimate normalization factors with calc_norm_factors.
4. Estimate dispersions and fit the model with differential_expression.
5. Convert the DEResult objects to a DataFrame or pass them to volcano and MA plotting helpers.

## Important implementation details
- _library_sizes and _row_totals summarize counts at the sample and gene level.
- _trimmed_weighted_mean is used for robust normalization.
- _solve_irls! updates the GLM coefficients during iterative fitting.
- The DataFrame and show methods make the results easy to inspect in notebooks or the REPL.

## Threading notes
- `calc_norm_factors()` now defaults to threaded sample-wise normalization when multiple threads are available.
- `estimate_dispersions()` now defaults to threaded gene-wise dispersion estimation.
- `combat_correction()` now defaults to threaded batch-wise correction, with each batch processed independently.

## Why this file matters
This module is the statistical backbone of transcriptomics analysis in BioToolkit. It turns raw count data into interpretable effect sizes and significance values that other modules can visualize or integrate into larger analyses.
