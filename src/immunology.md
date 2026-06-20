# `immunology.jl` - Immunology and Repertoire Analysis

## Overview

`immunology.jl` provides TCR/BCR repertoire parsing, clonotype grouping, VDJ/isotype summaries, diversity metrics, CDR3 physicochemical features, epitope-binding approximations, lineage trees, and HLA association helpers.

### Purpose

This page documents the public API in `immunology.jl`, grouped by the biological workflow it supports. Each entry describes the role of the function or type in practical BioToolkit analyses.

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

## 1. Repertoire Parsing

These functions normalize immune receptor records into clonotype- and segment-level summaries.

| API | Description |
|---|---|
| `extract_cdr3` | Extracts CDR3 sequences from receptor annotations or sequence records. |
| `clonotype_table` | Builds clonotype counts and frequencies. |
| `assign_vdj_segments` | Assigns V, D, and J segment labels from sequence/annotation evidence. |
| `paired_clonotype_table` | Combines paired-chain receptor information. |
| `igblast_assign` | Parses or approximates IgBLAST-style segment assignment. |
| `imgt_highvquest_manifest` | Builds an IMGT/HighV-QUEST-style manifest for submitted sequences. |

## 2. Diversity and Usage

Population-level summaries quantify repertoire diversity, expansion, and gene usage.

| API | Description |
|---|---|
| `isotype_switching_summary` | Summarizes B-cell isotype switching patterns. |
| `germline_usage_bias` | Tests or summarizes germline gene usage bias. |
| `repertoire_diversity_metrics` | Computes Shannon, Simpson, richness, clonality, and related metrics. |
| `clonal_expansion_test` | Tests for expanded clonotypes across groups. |
| `cdr3_length_spectrum` | Builds CDR3 length distributions. |
| `differential_vj_usage` | Tests differential V/J usage. |
| `v_gene_family_usage` | Aggregates V gene usage by family. |
| `somatic_hypermutation_rate` | Estimates BCR somatic hypermutation burden. |
| `bcr_affinity_maturation_score` | Scores affinity-maturation-like mutation patterns. |

## 3. Specificity, HLA, and Lineage

Advanced helpers score specificity, similarity, and lineage relationships.

| API | Description |
|---|---|
| `bepipred_like_scores` | Computes B-cell epitope propensity scores. |
| `mhcflurry_like_scores` | Scores peptide-MHC binding features. |
| `tcrdist_like_matrix` | Computes TCRdist-like pairwise receptor distances. |
| `antigen_specificity_clustering` | Clusters receptors by sequence/specificity features. |
| `clonotype_trajectory` | Tracks clonotype abundance over time or pseudotime. |
| `hla_type_enrichment` | Tests HLA enrichment across groups. |
| `cdr3_physicochemical_properties` | Computes CDR3 length, charge, hydrophobicity, and composition features. |
| `network_centrality_clonotypes` | Computes centrality metrics in clonotype similarity networks. |
| `convergent_clonotype_detection` | Detects similar clonotypes arising independently. |
| `lineage_tree_from_clones` | Builds lineage trees from related clone sequences. |
| `epitope_binding_profile` | Summarizes binding scores across epitopes/HLA alleles. |

---

## Complete Usage Example

```julia
using BioToolkit

clones = clonotype_table(receptor_records)
diversity = repertoire_diversity_metrics(clones)
expanded = clonal_expansion_test(clones, group_labels)
dist = tcrdist_like_matrix(clones)
```

