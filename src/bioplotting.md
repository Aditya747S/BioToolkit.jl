# bioplotting.jl

## Purpose
This file is the statistical plotting backbone of BioToolkit. It converts differential expression and GWAS-style analysis results into standardized plot records, then builds publication-ready figures from those records.

## Main structs
- VolcanoPoint: one point on a volcano plot, including fold change, p-value, adjusted p-value, category, and label flag.
- MAPoint: one point on an MA plot, storing abundance, fold change, p-value, adjusted p-value, and category.
- VolcanoPlotResult: the full volcano output with all points, a summary tuple, and the final figure.
- MAPlotResult: the full MA output with point records, summary data, and figure.
- ClusteredHeatmapResult: the reordered heatmap matrix together with row and column orderings and labels.
- ManhattanPoint: one GWAS point with chromosome, cumulative position, significance category, and label flag.
- ManhattanPlotResult: the Manhattan plot data and figure container.
- QQPoint: one QQ plot point with expected and observed -log10(p) values.
- QQPlotResult: the QQ plot point list and figure.
- ForestPoint: one meta-analysis or forest-plot row with effect size, standard error, confidence interval, and significance values.
- ForestPlotResult: the assembled forest plot output.

## Main public functions
- volcano_data(results; lfc_cutoff, fdr_cutoff, label_top): converts differential expression results into VolcanoPoint records and summary counts.
- volcano_plot(results; ...): builds the volcano figure and returns a VolcanoPlotResult.
- ma_data(results; lfc_cutoff, fdr_cutoff): converts differential expression results into MA-point records.
- ma_plot(results; ...): builds the MA figure and returns an MAPlotResult.
- clustered_heatmap(matrix; row_labels, column_labels, scale, metric): reorders the matrix by hierarchical clustering and returns a ClusteredHeatmapResult.
- manhattan_data(result; pvalue_threshold, label_top): converts GWAS results into ManhattanPoint records.
- manhattan_plot(result; ...): builds the Manhattan figure and returns a ManhattanPlotResult.
- qq_data(result): produces QQPoint records and genomic inflation summary statistics.
- qq_plot(result; ...): builds the QQ figure and returns a QQPlotResult.
- gwas_forest_plot(result): converts meta-analysis results into ForestPoint records and renders a forest plot.
- export_plot(result, save_path): writes the stored figure from a plot result to disk.

## What the module does
The module separates plotting into two phases. First, it normalizes raw analysis output into point records with consistent fields. Second, it builds a figure from those records using Plots.jl. This structure makes the plotting code reusable because the data preparation step can be inspected independently from the figure generation step.

## Typical workflows
For differential expression, a user usually calls volcano_plot or ma_plot on a collection of gene-level results. For GWAS, the user calls manhattan_plot, qq_plot, or gwas_forest_plot on a GWASResult or MetaAnalysisResult. For clustered heatmaps, the user passes a numeric matrix and optional row or column labels, and the module returns the reordered matrix together with the final figure.

## Important implementation details
- _safe_probability clamps invalid or missing-looking probabilities into a usable plotting range.
- _top_labels selects the strongest labels so plots do not become cluttered.
- _scale_matrix supports row-wise, column-wise, or no scaling.
- _hierarchical_order performs the clustering used to reorder heatmaps.
- _save_plots_figure handles optional file export.

## Threading notes
- `volcano_data()`, `ma_data()`, `manhattan_data()`, `qq_data()`, and `gwas_forest_plot()` now default to threaded point construction when multiple threads are available.
- `_scale_matrix()` now defaults to threaded row-wise or column-wise normalization.
- Figure construction remains serial; the threaded work is limited to independent data preparation.

## How to use it
1. Feed analysis results into the appropriate data function or plot function.
2. Use the returned summary to inspect overall significance and count statistics.
3. Reuse the stored figure directly or export it with export_plot.
4. If you need a custom plot, use the returned point records as a clean data source instead of re-running analysis.

## Why this file matters
This file keeps BioToolkit's statistical plots consistent across different domains. The point types and result types give the rest of the package a stable visualization interface, and the figures can be regenerated or restyled without recomputing the underlying biology.
