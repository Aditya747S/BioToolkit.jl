# `microbiome.jl` - Microbiome Community Analysis

## Overview

`microbiome.jl` provides community-profile containers, compositional transforms, alpha/beta diversity, UniFrac-style distances, ordination, differential abundance, co-occurrence networks, source tracking, taxonomic classification, MAG binning, strain profiling, and pathway summaries.

### Purpose

This page documents the public API in `microbiome.jl`, grouped by the biological workflow it supports. Each entry describes the role of the function or type in practical BioToolkit analyses.

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

Result types preserve ordination, differential abundance, network, and source-tracking metadata.

| API | Description |
|---|---|
| `CommunityProfile` | Taxa/features by samples abundance table with taxonomy/sample metadata. |
| `PCoAResult` | Principal coordinates result. |
| `NMDSResult` | NMDS ordination result. |
| `ANCOMResult` | ANCOM-style differential abundance result. |
| `SongbirdResult` | Songbird-style multinomial differential ranking result. |
| `MicrobiomeNetwork` | Co-occurrence network with taxa nodes and weighted edges. |
| `SourceTrackingResult` | Source contribution estimates for sink samples. |

## 2. Transforms and Diversity

Core functions compute compositional transforms and ecological distances/diversity.

| API | Description |
|---|---|
| `clr_transform` | Centered log-ratio transform. |
| `ilr_transform` | Isometric log-ratio transform. |
| `bray_curtis` | Bray-Curtis dissimilarity between samples. |
| `unifrac` | Unweighted UniFrac-style distance. |
| `weighted_unifrac` | Weighted UniFrac-style distance. |
| `pairwise_bray_curtis` | Pairwise Bray-Curtis distance matrix. |
| `pairwise_unifrac` | Pairwise UniFrac-style matrix. |
| `shannon_entropy` | Shannon alpha diversity. |
| `simpson_index` | Simpson diversity index. |
| `faith_pd` | Faith phylogenetic diversity. |

## 3. Ordination, Testing, and Networks

Higher-level functions support ordination plots, differential abundance, networks, and source tracking.

| API | Description |
|---|---|
| `pcoa` | Principal coordinates analysis. |
| `pcoa_plot` | Plot-ready PCoA payload. |
| `nmds` | Nonmetric multidimensional scaling. |
| `ancom` | ANCOM-style differential abundance. |
| `songbird` | Songbird-style differential ranking. |
| `cooccurrence_network` | Builds taxa co-occurrence network. |
| `network_plot` | Plot-ready network payload. |
| `source_tracking_model` | Fits a source-tracking model. |
| `source_tracking` | Estimates source contributions. |
| `source_tracking_posterior_summary` | Summarizes source-tracking posteriors. |
| `mag_bin_contigs` | Bins contigs into MAG-like groups. |
| `kraken_like_classify` | Classifies reads/contigs by k-mer votes. |
| `viral_contig_scores` | Scores contigs for viral signal. |
| `humann_like_pathways` | Aggregates taxa/gene families into pathway-like summaries. |
| `strainge_like_variants` | Profiles strain variants. |
| `lca_taxonomy_from_votes` | Computes lowest-common-ancestor taxonomy from votes. |
| `strain_haplotype_profile` | Builds strain haplotype profiles. |

---

## Complete Usage Example

```julia
using BioToolkit

profile = CommunityProfile(counts, taxa, samples)
dist = pairwise_bray_curtis(profile)
ord = pcoa(dist)
da = ancom(profile, group_labels)
net = cooccurrence_network(profile)
```

