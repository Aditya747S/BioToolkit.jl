# `spatial.jl` - Spatial Transcriptomics

## Overview

`spatial.jl` implements spatial experiment containers, reference matrix construction, RCTD/cell2location-like deconvolution, spatial variability tests, neighborhood graphs, domain clustering, ligand-receptor spatial scores, spot QC, coexpression modules, tissue boundary marking, spatial pseudotime, and several spatial-method approximations.

### Purpose

This documentation covers the public API implemented in `spatial.jl`, with concrete descriptions for the result types and workflow functions exposed by BioToolkit.

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

## 1. Types and Deconvolution

Spatial containers hold spot/cell coordinates, expression, metadata, and deconvolution outputs.

| API | Description |
|---|---|
| `SpatialExperiment` | Spatial expression container with coordinates and metadata. |
| `DeconvolutionResult` | Cell-type proportion/deconvolution result. |
| `build_reference_matrix` | Builds cell-type reference expression matrix. |
| `rctd_deconvolution` | Runs RCTD-like spot deconvolution. |
| `cell2location_deconvolution` | Runs cell2location-like deconvolution. |
| `cell2location_like_segmentation` | Segments spatial data from deconvolution-like outputs. |

## 2. Spatial Statistics and Domains

These functions detect spatially structured genes, neighborhoods, and domains.

| API | Description |
|---|---|
| `spatially_variable_genes` | Finds genes with spatial variability. |
| `spatial_autocorrelation` | Computes spatial autocorrelation statistics. |
| `spatial_neighborhood_graph` | Builds graph from spot/cell coordinates. |
| `spatial_domain_clustering` | Clusters tissue domains. |
| `bayesspace_like_domains` | BayesSpace-like spatial domain refinement. |
| `spagc_like_domains` | SpaGCN-like graph/domain clustering. |
| `sparkx_spatial_de` | SPARK-X-like spatial differential expression. |
| `spatialde_gp_de` | SpatialDE-like Gaussian-process spatial DE. |
| `spatial_markov_refine_domains` | Markov smoothing/refinement of spatial domains. |

## 3. Communication, QC, and Trajectories

Spatial downstream helpers score interactions, modules, boundaries, and trajectories.

| API | Description |
|---|---|
| `ligand_receptor_spatial_score` | Scores ligand-receptor spatial colocalization. |
| `niche_weighted_communication` | Weights communication by spatial niche/neighborhood. |
| `spatial_lr_permutation_test` | Permutation test for spatial ligand-receptor signal. |
| `spatial_trajectory_graph` | Builds trajectory graph over spatial neighborhoods. |
| `build_spot_deconvolution_qc` | Builds QC summaries for deconvolution outputs. |
| `spatial_coexpression_modules` | Finds spatial coexpression modules. |
| `mark_tissue_boundary_spots` | Flags tissue boundary spots. |
| `spatial_pseudotime` | Computes pseudotime over spatial graph. |

---

## Complete Usage Example

```julia
using BioToolkit

sp = SpatialExperiment(counts, coordinates)
ref = build_reference_matrix(singlecell_reference, labels)
dec = rctd_deconvolution(sp, ref)
svg = spatially_variable_genes(sp)
domains = spatial_domain_clustering(sp)
```

