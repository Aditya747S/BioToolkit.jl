# `singlecell_interop.jl` - Single-Cell Interoperability

## Overview

`singlecell_interop.jl` contains exchange helpers for moving BioToolkit single-cell objects to and from H5AD-style files and BioToolkit archive containers.

### Purpose

This page documents the public BioToolkit APIs implemented in `singlecell_interop.jl`. The sections follow the main workflow stages and describe how each function or type is intended to be used.

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

## 1. H5AD Exchange

H5AD helpers bridge BioToolkit experiments and AnnData-like files.

| API | Description |
|---|---|
| `read_h5ad` | Reads an H5AD-like file into a `SingleCellExperiment`. |
| `write_h5ad` | Writes a `SingleCellExperiment` to an H5AD-like file. |

## 2. Archive Exchange

Archive helpers store large experiments in reusable on-disk form.

| API | Description |
|---|---|
| `SingleCellArchive` | Archive container for assays, embeddings, layers, and metadata. |
| `save_singlecell_archive` | Writes a single-cell archive. |
| `load_singlecell_archive` | Loads archive metadata and handles. |
| `archive_to_singlecell_experiment` | Materializes an archive as a `SingleCellExperiment`. |

---

## Complete Usage Example

```julia
using BioToolkit

sce = read_h5ad("atlas.h5ad")
save_singlecell_archive(sce, "atlas.btarchive")
archive = load_singlecell_archive("atlas.btarchive")
sce2 = archive_to_singlecell_experiment(archive)
```

