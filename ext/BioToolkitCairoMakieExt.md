# BioToolkitCairoMakieExt.jl

## Purpose
This extension provides CairoMakie output for the genome-browser layer. It is the static-export companion to the generic Makie extension and is intended for file-based figure generation.

## Main responsibilities
- Build CairoMakie figures for `GenomeBrowser` objects.
- Dispatch `render!` for gene, coverage, and alignment tracks.
- Provide `export_figure(browser, path; dpi=300)` for direct file export.
- Keep export-only rendering separate from interactive display.

## Public behavior added by the extension
- `render!(axis, track::GeneTrack, viewport::GenomeViewport)` renders gene-track plans into CairoMakie axes.
- `render!(axis, track::CoverageTrack, viewport::GenomeViewport)` renders coverage-track plans.
- `render!(axis, track::AlignmentTrack, viewport::GenomeViewport)` renders alignment-track plans.
- `export_figure(browser, path; dpi=300)` saves a browser figure to disk and returns the output path.
- The figure layout mirrors the Makie extension so a browser looks the same regardless of whether it is displayed interactively or saved to disk.

## Internal helper
- `_browser_figure(browser)` constructs the CairoMakie figure from the browser viewport and tracks.

## How it is used
This extension is for users who want saved output rather than interactive display. The same browser render plans are reused, but CairoMakie handles the final raster/vector file output.

It is the preferred path for PDF, SVG, or PNG export in documentation and publication workflows.

## Why it matters
Genome-browser output often needs to be shared in papers, reports, or static documentation. This extension gives BioToolkit a lightweight export path without forcing CairoMakie into the base runtime.