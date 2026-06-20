# `singlecell.jl` - Single-Cell Analysis

## Overview

`singlecell.jl` provides SingleCellExperiment storage, normalization, SCTransform-like transforms, variable feature selection, PCA/UMAP/neighbors/clustering, markers, doublet detection, batch integration, spatial statistics, RNA velocity, ligand-receptor communication, WNN, perturbation prediction, H5AD/archive I/O, cell-type annotation, ambient RNA removal, and viewer helpers.

### Purpose

This page documents the public BioToolkit APIs implemented in `singlecell.jl`. The sections follow the main workflow stages and describe how each function or type is intended to be used.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Public API coverage** | Exported types and functions are grouped by use case. |
| **Concrete descriptions** | Entries describe real behavior rather than generic call placeholders. |
| **Analysis-ready outputs** | Results are suitable for downstream visualization, statistics, or reporting. |
| **Interoperable inputs** | APIs work with standard Julia matrices/tables and BioToolkit containers. |
| **Composable workflows** | Higher-level functions build on lower-level summaries where possible. |

---

## 1. Core Containers and Preprocessing

These APIs build and normalize single-cell experiment data.

| API | Description |
|---|---|
| `SingleCellExperiment` | Assays, row/cell metadata, reductions, graphs, and layers for single-cell data. |
| `SingleCellProjectionModel` | Model for projecting query cells into a reference space. |
| `count_matrix` | Retrieves count matrix from an experiment. |
| `normalize_counts` | Library-size/log normalizes counts. |
| `sctransform` | SCTransform-like variance-stabilizing normalization. |
| `find_variable_features` | Selects highly variable genes/features. |
| `gpu_singlecell_experiment` | Moves compatible arrays to GPU-backed storage. |
| `materialize_singlecell_experiment` | Converts lazy/GPU-backed data back to ordinary arrays. |

## 2. Embedding, Clustering, and Markers

Core analysis functions produce reductions, neighborhoods, clusters, pseudotime, and markers.

| API | Description |
|---|---|
| `fit_singlecell_projection_model` | Fits a reference projection model. |
| `project_singlecell` | Projects query cells into a model/reference. |
| `run_pca` | Computes PCA reductions. |
| `find_neighbors` | Builds nearest-neighbor graph. |
| `find_spatial_neighbors` | Builds spatial neighbor graph. |
| `run_umap` | Computes UMAP-like embedding. |
| `cluster_cells` | Clusters cells from graph or embedding. |
| `find_clusters` | Alias/helper for cluster discovery. |
| `calculate_pseudotime` | Computes pseudotime along a graph/trajectory. |
| `find_cluster_markers` | Finds markers per cluster. |
| `find_markers` | General marker testing. |
| `summarize_clusters` | Summarizes cluster sizes and metadata. |
| `cluster_marker_summary` | Combines cluster and marker summaries. |
| `detect_doublets` | Scores likely doublets. |
| `integrate_batches` | Integrates batches. |
| `integrate_data` | General integration wrapper. |
| `score_cell_cycle` | Scores cell-cycle phase. |

## 3. Spatial, Velocity, Communication, and Perturbation

Specialized single-cell workflows reuse experiment metadata, layers, embeddings, and graphs.

| API | Description |
|---|---|
| `SpatialMoranResult` | Spatial autocorrelation result. |
| `attach_spatial_coords!` | Adds spatial coordinates to cells/spots. |
| `moran_i_test` | Tests Moran I spatial autocorrelation. |
| `find_spatially_variable_genes` | Finds genes with spatial structure. |
| `RNAVelocityResult` | RNA velocity result. |
| `DynamicalRNAVelocityResult` | Dynamical velocity model result. |
| `attach_velocity_layers!` | Adds spliced/unspliced or velocity layers. |
| `calculate_rna_velocity` | Computes RNA velocity vectors. |
| `calculate_dynamical_rna_velocity` | Fits dynamical RNA velocity summaries. |
| `LigandReceptorPair` | Ligand-receptor database entry. |
| `CellCommunicationResult` | Ligand-receptor test result. |
| `CellCommunicationNetwork` | Cell-type communication graph. |
| `CommunicationPathwaySummary` | Pathway-level communication summary. |
| `LigandReceptorReport` | Ranked LR report. |
| `default_ligand_receptor_pairs` | Built-in LR pair list. |
| `find_cell_communication` | Tests ligand-receptor communication. |
| `communication_network` | Builds communication network. |
| `communication_pathway_summary` | Summarizes pathway-level communication. |
| `rank_ligand_receptor_report` | Ranks LR pairs. |
| `WNNResult` | Weighted nearest-neighbor result. |
| `weighted_nearest_neighbors` | Combines modalities into WNN graph. |
| `PerturbationPredictionResult` | Perturbation prediction result. |
| `predict_perturbation` | Predicts perturbation effects. |
| `CellTypeAnnotationResult` | Cell type annotation result. |
| `annotate_cell_types` | Annotates cells from markers/reference. |
| `AmbientRNARemovalResult` | Ambient RNA correction result. |
| `remove_ambient_rna` | Removes ambient RNA contamination. |

## 4. I/O, Archives, Plots, and Viewer

Storage and visualization helpers keep large experiments portable.

| API | Description |
|---|---|
| `save_singlecell_experiment` | Saves an experiment. |
| `load_singlecell_experiment` | Loads an experiment. |
| `SingleCellArchive` | On-disk archive container. |
| `save_singlecell_archive` | Writes archive data. |
| `load_singlecell_archive` | Loads archive metadata/layers. |
| `archive_to_singlecell_experiment` | Materializes an archive as an experiment. |
| `read_h5ad` | Reads AnnData/H5AD-like files. |
| `write_h5ad` | Writes H5AD-like files. |
| `plot_trajectory` | Plot-ready trajectory payload. |
| `plot_rna_velocity` | Plot-ready velocity payload. |
| `plot_rna_velocity_quiver` | Quiver payload for velocity vectors. |
| `plot_communication_network` | Communication network plot payload. |
| `plot_communication_pathways` | Pathway plot payload. |
| `plot_ligand_receptor_report` | LR report plot payload. |
| `SingleCellViewer` | Interactive viewer state. |
| `interactive_singlecell_viewer` | Creates viewer from an experiment. |
| `lasso_select_cells` | Selects cells in embedding space. |
| `recluster_singlecell_viewer!` | Reclusters viewer selection. |
| `cell_hover_text` | Formats hover labels. |
| `top_expressed_genes` | Returns top genes for cells/clusters. |

---

## Complete Usage Example

```julia
using BioToolkit

sce = SingleCellExperiment(counts)
sce = normalize_counts(sce)
sce = find_variable_features(sce)
sce = run_pca(sce)
sce = run_umap(sce)
clusters = cluster_cells(sce)
```

