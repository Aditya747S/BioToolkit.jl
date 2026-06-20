# `enrichment.jl` - Gene Set Enrichment

## Overview

`enrichment.jl` manages annotation databases, identifier mapping, over-representation analysis, GSEA/GSVA-style scoring, pathway collections, network propagation, and plot-ready enrichment summaries.

### Purpose

This file documents the exported analysis objects and workflow functions in `enrichment.jl`. It focuses on what each API does, what stage of the workflow it belongs to, and how the pieces compose with the rest of BioToolkit.

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

## 1. Databases and IDs

Annotation databases map terms to genes and translate identifiers across namespaces.

| API | Description |
|---|---|
| `IDMapper` | Identifier mapping table/container. |
| `EnrichmentTerm` | One gene-set term with id, name, namespace, and genes. |
| `EnrichmentDatabase` | Collection of enrichment terms and optional ID mappings. |
| `EnrichmentResult` | Over-representation result for one term. |
| `GSEAResult` | Rank-based gene-set enrichment result. |
| `GSVAResult` | Sample-level gene-set activity result. |
| `load_annotation_database` | Loads a serialized annotation database. |
| `save_annotation_database` | Saves an annotation database. |
| `build_annotation_database` | Builds a database from term-to-gene records. |
| `builtin_annotation_database` | Loads built-in term collections. |
| `builtin_annotation_terms` | Lists built-in terms. |
| `map_id` | Maps one identifier. |
| `map_ids` | Maps many identifiers. |

## 2. Enrichment Tests

ORA, GSEA, GSVA, and collection-specific wrappers share result formats.

| API | Description |
|---|---|
| `enrichment_test` | Runs hypergeometric/Fisher-style over-representation analysis. |
| `go_enrichment` | Runs GO enrichment using a GO-like database. |
| `kegg_enrichment` | Runs KEGG enrichment. |
| `fgsea_like` | Fast preranked GSEA-like scoring. |
| `gsea_preranked` | Runs rank-based enrichment from ordered genes/statistics. |
| `gsva_score` | Computes sample-level pathway activity scores. |
| `reactome_enrichment` | Runs Reactome-style enrichment. |
| `wikipathways_enrichment` | Runs WikiPathways-style enrichment. |
| `msigdb_enrichment` | Runs MSigDB-style enrichment. |
| `competitive_gene_set_test` | Tests gene sets against background genes. |
| `self_contained_gene_set_test` | Tests whether a gene set has signal without comparing to other genes. |

## 3. Networks and Reporting

Reporting helpers build overlaps, maps, semantic summaries, and plot payloads.

| API | Description |
|---|---|
| `dotplot` | Creates dotplot-ready enrichment data. |
| `network_propagation` | Propagates gene scores across a network. |
| `heat_diffusion_enrichment` | Scores terms after heat-diffusion smoothing. |
| `gene_set_overlap_matrix` | Computes pairwise term overlaps. |
| `jaccard_similarity_matrix` | Computes pairwise Jaccard similarities. |
| `leading_edge_genes` | Extracts leading-edge genes from GSEA results. |
| `enrichment_map` | Builds an enrichment-map graph. |
| `enrichment_heatmap_data` | Builds heatmap-ready gene-set matrices. |
| `bubble_chart_data` | Builds bubble-chart payloads. |
| `gene_ontology_semantic_similarity` | Computes GO semantic similarity between terms. |
| `go_slim_mapping` | Maps detailed GO terms to slim terms. |
| `term_to_gene_matrix` | Builds a binary term-by-gene matrix. |
| `rank_genes_by_set_membership` | Ranks genes by membership frequency or weighted set scores. |

---

## Complete Usage Example

```julia
using BioToolkit

db = builtin_annotation_database(:go)
res = go_enrichment(significant_genes, db; universe=background_genes)
plotdata = dotplot(res)
gsea = gsea_preranked(ranked_gene_scores, db)
```

