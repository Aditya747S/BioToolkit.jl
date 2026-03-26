# singlecell.jl

## Purpose
This file implements a single-cell RNA-seq analysis workflow. It covers count storage, normalization, variance stabilization, dimensionality reduction, clustering, marker detection, doublet detection, batch correction, and cell-cycle scoring.

## Main type
- `SingleCellExperiment` stores sparse counts, gene and cell IDs, metadata, dimensionality reductions, and cluster assignments.

## Public functions
- Input and preprocessing: `count_matrix`, `normalize_counts`, and `sctransform`.
- Dimensionality reduction: `run_pca` and `run_umap`.
- Clustering and summaries: `cluster_cells`, `summarize_clusters`, `cluster_marker_summary`, and `find_cluster_markers`.
- Quality/control workflows: `detect_doublets`, `integrate_batches`, and `score_cell_cycle`.

## How it is used
Users typically build a `SingleCellExperiment` from a matrix plus gene and cell labels, normalize the counts, run `run_pca`, and then cluster cells with `cluster_cells`. Marker detection uses `find_cluster_markers` or `cluster_marker_summary`, while `detect_doublets` and `integrate_batches` help clean up and harmonize samples before final interpretation.

`score_cell_cycle` gives a lightweight phase label based on predefined S and G2M gene sets, and `run_umap` provides a visual embedding either through UMAP when available or a fallback two-dimensional projection.

## Implementation notes
- The object stores intermediate analysis artifacts directly on the experiment instance so downstream steps can reuse them.
- `cluster_cells` supports KNN, graph-based, and k-means-style clustering modes.
- `run_umap` will use a UMAP implementation in `Main` when present and otherwise falls back gracefully.
- `detect_doublets` simulates synthetic doublets in the embedding space and scores cells by nearest-neighbor composition.

## Why it matters
Single-cell analysis involves a long sequence of preprocessing and interpretation steps. This file keeps that sequence together in one reusable container-based API, which makes it easier to build consistent workflows across datasets.
