# `gwas.jl` - GWAS, Genetics, and Fine Mapping

## Overview

`gwas.jl` implements genotype I/O, association scans, sample/variant QC, LD and GRM utilities, PRS, meta-analysis, fine mapping, Mendelian randomization, heritability, colocalization, ancestry, liftover, simulation, and result persistence.

### Purpose

This document is a source-matched reference for the public API in `gwas.jl`. It groups exported types and functions by workflow stage and explains the biological or statistical role of each entry.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Workflow-oriented grouping** | Functions are organized by how users run the analysis. |
| **Typed domain state** | Reusable records and result objects are represented explicitly. |
| **Scalable summaries** | APIs return compact matrices/tables for large biological datasets. |
| **Interoperable outputs** | Results can feed plotting, enrichment, reporting, and downstream modeling modules. |
| **Explicit thresholds and controls** | Filtering, QC, and model functions expose important parameters as keywords. |

---

## 1. Types and I/O

Core containers and readers/writers handle PLINK, BGEN, and summary-statistic style data.

| API | Description |
|---|---|
| `GenotypeMatrix` | Samples-by-variants genotype/dosage matrix with variant and sample metadata. |
| `GWASResult` | Association result table plus metadata. |
| `MetaAnalysisResult` | Combined-effect meta-analysis result. |
| `BedReader` | PLINK BED reader state. |
| `BgenReader` | BGEN reader state. |
| `BgenVariant` | One BGEN variant metadata record. |
| `read_plink` | Reads PLINK BED/BIM/FAM files. |
| `write_plink` | Writes PLINK-style genotype files. |
| `read_bed_genotypes` | Reads BED genotype matrix data. |
| `read_bgen` | Reads BGEN variants/dosages. |
| `write_bgen` | Writes BGEN-like dosage data. |
| `to_plink_dataframe` | Converts genotype/variant data to PLINK-like tables. |
| `write_plink_sumstats` | Writes summary statistics in PLINK-compatible form. |
| `save_gwas_result` | Serializes a GWAS result. |
| `load_gwas_result` | Loads a serialized GWAS result. |

## 2. Association Testing

Scans support common phenotypes and model families.

| API | Description |
|---|---|
| `gwas_linear_scan` | Linear regression GWAS for quantitative traits. |
| `gwas_lmm_scan` | Linear mixed-model scan with relatedness correction. |
| `gwas_logistic_scan` | Logistic regression GWAS for binary traits. |
| `gwas_gxe_interaction` | Gene-by-environment interaction scan. |
| `score_test_linear` | Score test for linear models. |
| `likelihood_ratio_test` | Likelihood-ratio model comparison. |
| `gwas_survival_scan` | Survival phenotype association scan. |
| `gwas_ordinal_scan` | Ordinal phenotype scan. |
| `gwas_multivariate_scan` | Multivariate phenotype scan. |
| `burden_test` | Region/gene burden test. |
| `skat_test` | SKAT-style variance-component test. |
| `permutation_pvalue` | Empirical p-value by permutation. |
| `gwas_power_calculation` | Power/sample-size utility. |

## 3. QC, LD, Relatedness, and Ancestry

QC functions operate on variants, samples, LD, and relatedness structure.

| API | Description |
|---|---|
| `calculate_maf` | Computes minor allele frequencies. |
| `calculate_missingness` | Computes variant/sample missingness. |
| `missingness_filter` | Filters variants/samples by missingness. |
| `info_score_filter` | Filters imputed variants by INFO score. |
| `hwe_exact` | Exact Hardy-Weinberg test. |
| `calculate_hwe_pvalues` | Computes HWE p-values across variants. |
| `hwe_filter` | Filters variants by HWE p-value. |
| `genomic_control_lambda` | Computes genomic inflation factor. |
| `apply_genomic_control` | Adjusts statistics by genomic control. |
| `calculate_ld_matrix` | Computes pairwise LD. |
| `calculate_grm` | Computes genomic relationship matrix. |
| `compute_ld_scores` | Computes LD scores. |
| `prune_ld` | LD pruning. |
| `ld_clumping` | P-value guided LD clumping. |
| `calculate_ibd` | IBD estimate. |
| `calculate_king_kinship` | KING-like kinship. |
| `detect_related_pairs` | Finds related sample pairs. |
| `sample_prune_relatedness` | Prunes samples by relatedness. |
| `ancestry_inference` | Infers ancestry from genotype PCs/reference data. |
| `gwas_pca` | Fits genotype PCA. |
| `project_pca` | Projects samples onto PCA loadings. |
| `calculate_sex_check` | Sex-check summary from genotype data. |
| `calculate_heterozygosity_outliers` | Finds heterozygosity outliers. |
| `sample_qc_report` | Builds a sample QC report. |
| `filter_variants` | General variant filtering. |
| `filter_samples` | General sample filtering. |
| `merge_genotype_matrices` | Merges compatible genotype matrices. |
| `flip_alleles` | Flips allele coding. |
| `harmonise_alleles` | Harmonizes alleles between datasets. |
| `normalise_chromosome` | Normalizes chromosome labels. |

## 4. PRS, Meta-Analysis, Fine Mapping, and Causal Genetics

Advanced workflows reuse association results and LD reference data.

| API | Description |
|---|---|
| `calculate_prs` | Computes polygenic risk scores. |
| `prs_ldpred` | LDpred-like PRS shrinkage. |
| `prs_cross_validation` | Cross-validates PRS settings. |
| `meta_analyze` | Fixed/random-effects meta-analysis. |
| `overlap_gwas_peaks` | Finds overlapping significant loci. |
| `gene_based_test` | Aggregates variant signal to genes. |
| `LDSCResult` | LD score regression result. |
| `ldsc_heritability` | SNP-heritability by LDSC. |
| `estimate_heritability_greml` | GREML-style heritability estimate. |
| `partitioned_heritability` | Partitioned heritability by annotation. |
| `ldsc_genetic_correlation` | Genetic correlation by LDSC. |
| `SuSiEResult` | SuSiE fine-mapping result. |
| `fine_map_susie` | Runs SuSiE-style fine mapping. |
| `calculate_credible_set` | Builds credible sets from posterior probabilities. |
| `posterior_inclusion_probability` | Computes PIP values. |
| `MRResult` | Mendelian randomization result. |
| `mr_two_sample` | Two-sample MR estimate. |
| `mr_egger` | MR-Egger regression. |
| `mr_pleiotropy_test` | Tests directional pleiotropy. |
| `conditional_analysis` | Conditioned association analysis. |
| `cojo_stepwise` | COJO-style stepwise conditional selection. |
| `joint_analysis` | Joint multi-variant model. |
| `coloc_result` | Colocalization result container. |
| `coloc_abf` | Approximate Bayes factor colocalization. |
| `conditional_fdr` | Conditional FDR analysis. |
| `storey_pi0_estimate` | Storey pi0 estimate for q-values. |
| `genomic_sem_fit` | Genomic SEM-like model fit. |
| `twas_scan` | Transcriptome-wide association scan. |
| `smr_test` | SMR-style test. |
| `popcorn_genetic_correlation` | Cross-ancestry genetic correlation estimate. |
| `ebi_lookup` | EBI/GWAS lookup helper. |
| `liftover` | Coordinate liftover helper. |
| `functional_annotation` | Annotates variants with functional categories. |
| `ld_expand_credible_set` | Expands credible sets by LD. |
| `estimate_effective_n` | Effective sample size estimate. |
| `simulate_genotypes` | Simulates genotype matrices. |
| `simulate_phenotype` | Simulates phenotypes. |
| `dosage_to_hardcall` | Converts dosages to hard calls. |
| `info_score_from_dosage` | Computes imputation INFO from dosages. |
| `rank_inverse_normal` | Rank inverse-normal transform. |
| `rank_inverse_normal!` | In-place rank inverse-normal transform. |

---

## Complete Usage Example

```julia
using BioToolkit

gt = read_plink("cohort")
qc = gwas_qc_report(gt)
res = gwas_linear_scan(gt, phenotype)
res = apply_genomic_control(res)
prs = calculate_prs(gt, res)
```

