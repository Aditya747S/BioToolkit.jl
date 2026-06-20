# `browser.jl` - Interactive Browser Utilities

## Overview

`browser.jl` contains lightweight in-memory browser/viewer state for single-cell archives and interactive selection workflows.

### Purpose

This page documents the public surface of `browser.jl`: the result types, workflow functions, and helper APIs a BioToolkit user is expected to call directly. Private helpers are intentionally omitted unless they define behavior users must understand.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Structured domain outputs** | Stable results are returned as structs, named tuples, or table-friendly payloads. |
| **Composable Julia inputs** | APIs favor ordinary arrays, matrices, dictionaries, BioToolkit records, and DataFrames. |
| **Workflow granularity** | Each public function maps to a recognizable analysis or data-preparation step. |
| **Optional external integration** | Wrapper-style functions prepare command inputs or parse results without making all workflows depend on external tools. |
| **Reporting-ready summaries** | Many functions return data that can be plotted, exported, or passed into downstream modules. |

---

## 1. Archive Browsing

Archive helpers expose layers, embeddings, and metadata without forcing users to materialize every matrix.

| API | Description |
|---|---|
| `SingleCellArchiveBrowser` | State object for browsing a `SingleCellArchive`. |
| `browse_singlecell_archive` | Creates a browser object from an archive path or loaded archive. |
| `archive_layer` | Retrieves a named matrix layer from an archive. |
| `archive_embedding` | Retrieves a named low-dimensional embedding from an archive. |

## 2. Interactive Single-Cell Viewer

Viewer functions keep selection and reclustering state in Julia objects.

| API | Description |
|---|---|
| `SingleCellViewer` | Interactive single-cell viewer state with embedding, labels, selections, and metadata. |
| `interactive_singlecell_viewer` | Constructs a viewer from a `SingleCellExperiment`. |
| `lasso_select_cells` | Selects cells inside a polygon/lasso region in embedding space. |
| `recluster_singlecell_viewer!` | Reclusters the current viewer selection or full viewer dataset in place. |
| `cell_hover_text` | Formats metadata-rich hover labels for cells. |
| `top_expressed_genes` | Returns top marker/expression genes for selected cells or clusters. |

---

## Complete Usage Example

```julia
using BioToolkit

viewer = interactive_singlecell_viewer(sce; embedding=:umap)
selected = lasso_select_cells(viewer, polygon_points)
labels = cell_hover_text(viewer, selected)
```

