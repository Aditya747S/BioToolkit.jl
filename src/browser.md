# browser.jl

## Purpose
This file implements the genome browser core for BioToolkit. It defines the browser viewport, track types, render-plan objects, and the logic that turns genomic content into layouts that can be rendered by different graphical backends.

## Main structs
- GenomeViewport: stores the chromosome, genomic range, pixel width, and derived resolution.
- AbstractTrack: abstract parent for all browser tracks.
- GeneTrack: gene annotation track with style, height, color, and labeling settings.
- CoverageTrack: coverage or signal track with windowing and sorting settings.
- AlignmentTrack: read-alignment track with mismatch and pair display options.
- GenomeBrowser: the full browser object holding tracks, viewport, spacing, and title.
- GeneSegment: one drawn segment inside a gene placement.
- GenePlacement: a gene interval assigned to a row, with label and segment breakdown.
- GeneRenderPlan: complete gene-track layout for a viewport.
- CoverageBin: one binned coverage interval with summary statistics.
- CoverageRenderPlan: complete coverage-track layout for a viewport.
- ReadMismatch: one mismatch position inside a read.
- ReadPlacement: a read assigned to a row with blocks and mismatch annotations.
- ReadConnector: a paired-read connector line.
- AlignmentRenderPlan: complete alignment-track layout for a viewport.

## Public functions and constructors
- GenomeViewport(chrom, range, pixel_width): create a viewport from a genomic interval.
- GenomeBrowser(tracks, viewport; gap, title): assemble a browser from tracks and an existing viewport.
- GenomeBrowser(tracks, chrom, start, stop; width, gap, title): convenience constructor that builds the viewport for you.
- GeneTrack(intervals; style, height, color, label_field, max_labels): create a gene track from genomic intervals.
- GeneTrack(collection; kwargs...): construct a gene track from an interval collection.
- CoverageTrack(source; chrom, height, color, style, window_size, sorted, start): configure a coverage track.
- AlignmentTrack(source; show_mismatches, height, color, max_reads, rng_seed, show_pairs): configure a read-alignment track.
- genome_lod(viewport): decide whether the viewport should be shown at chromosome, gene, or basepair detail.
- render_track(track, viewport): convert one track into its corresponding render plan.
- render_browser(browser): convert all browser tracks into render plans.
- export_figure(browser, path; dpi): export browser figures when a graphics backend is available.

## What the module does
The browser code does not directly draw shapes. Instead, it computes an appropriate representation of what should be drawn. For example, genes can be packed into rows, coverage can be binned, and alignments can be downsampled or expanded depending on the current resolution. The result is a structured render plan that other code can turn into a figure.

## Rendering logic
The module uses level-of-detail thresholds to decide how much information to display. At coarse resolution, tracks are summarized. At intermediate resolution, genes are laid out by row. At fine resolution, alignments and read mismatches are shown in more detail. This keeps the browser readable across a wide range of genomic spans.

## Typical usage
1. Create a GenomeViewport for the chromosome and region of interest.
2. Create one or more tracks, such as GeneTrack, CoverageTrack, or AlignmentTrack.
3. Build a GenomeBrowser from the tracks and viewport.
4. Call render_browser to inspect the plans or use the plotting extension to render the browser visually.

## Important internal helpers
- _visible_interval and _visible_intervals filter features to the viewport.
- _gene_label chooses a readable label from feature metadata.
- _coerce_segment and _gene_segments break feature locations into drawable pieces.
- _pack_genes assigns genes to rows so overlapping annotations do not collide.
- _coverage_plan and _alignment_plan build the respective render plans.

## Why this file matters
This file is the architectural center of BioToolkit's genome browser. It keeps the layout rules in one place, so the same browser object can be rendered in different front ends without duplicating the logic for packing genes or summarizing reads.
