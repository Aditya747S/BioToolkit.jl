# BioToolkitPlotsExt.jl

## Purpose
This extension connects BioToolkit result objects to Plots.jl. It provides concrete `plot` methods for statistical, clinical, phylogenetic, and genome-browser outputs so callers can render figures without wiring each result type by hand.

## Main responsibilities
- Render BioToolkit differential-expression, GWAS, and clinical result objects with Plots.jl.
- Provide a Plots-based browser renderer for gene, coverage, and alignment tracks.
- Provide a Plots-based tree renderer for `PhyloTree`.
- Add direct `plot(...)` dispatch for BioToolkit result containers so existing plotting workflows keep working.
- Keep the plotting layer backend-specific while leaving result construction inside the core modules.

## Key helpers
- `_save_plot(path, figure)` writes a figure when an output path is supplied.
- `_km_step_data(result)` converts Kaplan-Meier survival output into stepwise x/y coordinates.
- `_browser_track_plot(plan, track; title="")` renders gene-browser plans.
- `_roc_figure(result; title="Survival ROC")` renders ROC curves with the no-skill baseline.
- `_cif_figure(result; title="Competing risks cumulative incidence")` renders competing-risk cumulative incidence curves.

## Public behavior added by the extension
- `volcano_plot`, `ma_plot`, `clustered_heatmap`, `manhattan_plot`, `qq_plot`, and `gwas_forest_plot` build Plots-backed figures from BioToolkit analysis results.
- `kaplan_meier_plot` renders one or many Kaplan-Meier curves and can overlay censor markers.
- `Plots.plot(result::VolcanoPlotResult)` and similar methods return the stored figure for already-built plot result objects.
- `Plots.plot(result::ROCResult)`, `Plots.plot(result::CIFResult)`, `Plots.plot(result::CoxResult)`, and `Plots.plot(result::OncoprintResult)` provide clinical plotting shortcuts.
- `Plots.plot(browser::GenomeBrowser)` renders stacked browser tracks using the browser render plans.
- `Plots.plot(plan::GeneRenderPlan)`, `Plots.plot(plan::CoverageRenderPlan)`, and `Plots.plot(plan::AlignmentRenderPlan)` expose direct plan rendering.
- `Plots.plot(tree::PhyloTree)` draws a rooted tree from `coordinates(tree)` and terminal labels.
- `pair_plot`, `violin_plot`, and `mosaic_plot` provide additional matrix-style and grouped-data visualizations.
- `Plots.plot(result::KaplanMeierResult)` and `Plots.plot(results::AbstractVector{<:KaplanMeierResult})` provide survival-curve shortcuts.
- `Plots.plot(result::ManhattanPlotResult)`, `Plots.plot(result::QQPlotResult)`, and `Plots.plot(result::ForestPlotResult)` expose the stored GWAS figures directly.

## Plot design notes
- Volcano and MA plots use cutoff lines and category-specific palettes so significance is visible at a glance.
- Manhattan plots sort chromosomes into a stable genome order and alternate chromosome colors.
- QQ plots draw a 45-degree reference line plus a beta-distribution confidence envelope.
- Forest plots sort hits by q-value and draw confidence intervals as horizontal error bars.
- Genome-browser plots preserve the render-plan structure from the core browser module so the visuals match the browser logic rather than inventing new layout rules.
- Tree plots rely on the phylogeny coordinates map to keep branch lengths and leaf order consistent.

## How it is used
This extension is activated when Plots.jl is available. Users can call `plot(result)` on a BioToolkit result object or browser/tree container and receive a ready-to-display figure.

The extension is deliberately thin: it keeps the styling and layout logic close to the plotting backend while leaving the scientific data preparation inside the core modules.
It is also the place where export-friendly plot objects are converted into a displayable `Figure` without duplicating the statistics that produced the result.

## Why it matters
BioToolkit's core modules produce rich analysis objects. This extension turns those objects into usable figures without forcing each analysis layer to duplicate plotting code.