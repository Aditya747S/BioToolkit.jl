# gwas.jl

## Purpose
This file implements the GWAS analysis layer of BioToolkit. It covers genotype storage, PLINK file I/O, association scans, polygenic scoring, LD clumping, meta-analysis, and downstream overlap and gene-based summaries.

## Main structs
- GenotypeMatrix: a decoded genotype matrix together with raw PLINK bytes, BIM metadata, FAM metadata, and an optional file prefix.
- GWASResult: association results for one study or phenotype, including SNP IDs, chromosomes, positions, alleles, effect sizes, standard errors, Z-scores, p-values, sample size, covariates, phenotype name, and method name.
- MetaAnalysisResult: multi-study summary statistics with q-values and heterogeneity metrics such as tau² and I².

## Public constructors and functions
- GenotypeMatrix(decoded, bim, fam; prefix): build a genotype matrix from decoded values.
- GWASResult(...): normalize GWAS result vectors into pooled, typed fields.
- MetaAnalysisResult(...): normalize meta-analysis outputs into a compact typed result.
- read_plink(prefix): load .bed, .bim, and .fam files into a GenotypeMatrix.
- write_plink(prefix, genotypes, bim, fam): write PLINK-compatible files back to disk.
- gwas_linear_scan: run a linear association test.
- gwas_lmm_scan: run a mixed-model association test.
- calculate_prs, prs_ldpred, prs_cross_validation: polygenic risk score routines.
- ld_clumping: reduce correlated hits by linkage disequilibrium.
- meta_analyze: combine association results across cohorts.
- overlap_gwas_peaks: intersect GWAS hits with peak sets.
- gene_based_test: aggregate evidence by gene.

## What the module does
The module takes users from raw genotype files to statistical association results. It can read PLINK data into an in-memory matrix, perform variant-level scans, filter and summarize hits, and then combine studies when more than one cohort is available. The result objects are designed to be compact and easy to plot or pass into downstream interpretation layers.

## How the structs work together
GenotypeMatrix is the input container that ties together genotypes and file metadata. GWASResult stores one scan's worth of statistics in a form that can be plotted, exported, or intersected with genomic regions. MetaAnalysisResult expands that same idea to combined results from multiple studies.

## Typical usage
1. Read a PLINK prefix with read_plink.
2. Run gwas_linear_scan or gwas_lmm_scan to obtain a GWASResult.
3. Apply ld_clumping or overlap_gwas_peaks to reduce or interpret the hits.
4. Use calculate_prs or prs_ldpred when you need risk-score style outputs.
5. Call meta_analyze to combine multiple studies into a MetaAnalysisResult.

## Important implementation details
- _read_bim and _read_fam parse the metadata tables into DataFrames.
- _decode_bed and _encode_bed translate PLINK binary genotypes to and from dense matrices.
- _covariate_matrix builds regression design matrices.
- _project_out and other helpers support model fitting and residualization.
- The result constructors convert inputs into pooled arrays and typed vectors to keep memory usage low.

## Why this file matters
This is the main GWAS analysis surface in BioToolkit. It connects genotype storage, statistical testing, and downstream interpretation in a single module.
