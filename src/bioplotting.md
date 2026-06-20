# `bioplotting.jl` - Plot-Ready Bioinformatics Data

## Overview

`bioplotting.jl` converts analysis results into stable, plot-ready data structures for volcano plots, MA plots, clustered heatmaps, Manhattan plots, QQ plots, and forest plots.

### Purpose

This page documents the public surface of `bioplotting.jl`: the result types, workflow functions, and helper APIs a BioToolkit user is expected to call directly. Private helpers are intentionally omitted unless they define behavior users must understand.

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

## 1. Result Types

Typed containers keep plotted points and plot metadata together.

| API | Description |
|---|---|
| `VolcanoPoint` | One volcano-plot point with effect size, significance, label, and optional grouping. |
| `VolcanoPlotResult` | Collection of volcano points plus thresholds and labels. |
| `MAPoint` | One MA-plot point with mean abundance and log-ratio. |
| `MAPlotResult` | Collection of MA points and plot metadata. |
| `ClusteredHeatmapResult` | Heatmap matrix with row/column ordering and clustering metadata. |
| `ManhattanPoint` | One genomic association point with chromosome, position, p-value, and label. |
| `ManhattanPlotResult` | Collection of Manhattan points with chromosome layout metadata. |
| `QQPoint` | Observed/expected p-value point for QQ plots. |
| `QQPlotResult` | QQ points plus genomic-control and confidence-band metadata. |
| `ForestPoint` | One forest-plot estimate with interval and label. |
| `ForestPlotResult` | Collection of forest points for meta-analysis or Cox/GWAS summaries. |

## 2. Plot Builders

Builders accept domain results and return data; rendering can be handled by extensions.

| API | Description |
|---|---|
| `volcano_data` | Creates `VolcanoPoint`s from differential-expression or generic result tables. |
| `volcano_plot` | Builds a complete volcano plot payload. |
| `ma_data` | Creates MA-plot points from expression/abundance results. |
| `ma_plot` | Builds a complete MA plot payload. |
| `clustered_heatmap` | Clusters rows/columns and returns ordered matrix payloads. |
| `manhattan_data` | Converts GWAS result tables into Manhattan points. |
| `manhattan_plot` | Builds a full Manhattan plot payload. |
| `qq_data` | Computes expected and observed p-values for QQ diagnostics. |
| `qq_plot` | Builds a full QQ plot payload. |
| `gwas_forest_plot` | Builds forest-plot intervals for loci, cohorts, or meta-analysis rows. |
| `export_plot` | Writes plot payloads through available plotting/export backends. |

---

## Complete Usage Example

```julia
using BioToolkit

volcano = volcano_data(de_results; logfc=:log2foldchange, pvalue=:padj)
qq = qq_plot(gwas_results)
forest = gwas_forest_plot(meta_results)
```

