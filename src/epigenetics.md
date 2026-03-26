# epigenetics.jl

## Purpose
This file implements the epigenomics analysis layer of BioToolkit. It covers chromatin accessibility, peak support, methylation analysis, and 3D genome structure objects, including single-cell chromatin support and contact-based topology features.

## Main structs
- SparseCoverageVector: sparse depth information for a genomic span.
- Peak: one peak record with coordinates, summit, score, and significance values.
- PeakSet: a collection of peaks plus chromosome lookup tables.
- PeakSupport: per-peak support metrics across samples.
- Epigenome: a combined epigenomic container with intervals, coverage, counts, and sample metadata.
- MethylationCall: one methylation observation at a genomic position in a sample.
- MethylationExperiment: methylated and total count matrices with regions and sample metadata.
- MethylationResult: one differential methylation result row.
- SingleCellChromatinExperiment: single-cell chromatin container with open peaks, reductions, and metadata.
- CoaccessibilityEdge: one edge between two peaks with a correlation score.
- TadResult: one topologically associating domain result with genomic bin bounds and state.
- ContactMatrix: a sparse contact matrix with genomic bins and chromosome annotation.

## Public functions
- coverage_depth(coverage, position): query depth at one genomic position.
- coverage_segments(coverage; chrom): convert sparse coverage into coverage segments.
- calculate_coverage, coverage_depth, coverage_segments: coverage-related helpers.
- call_peaks: identify peak regions from chromatin signal.
- normalize_gc_bias: adjust signal for GC-related bias.
- count_overlaps and summarize_peak_support: support peak-by-sample summarization.
- differential_binding: compare peak accessibility between groups.
- bin_methylation and differential_methylation: methylation analysis helpers.
- tfidf, run_lsi, rsvd, gene_activity_score: single-cell chromatin dimensionality and feature transforms.
- calculate_coaccessibility: estimate peak co-accessibility.
- compute_motif_deviations and detect_footprints: motif and footprinting analysis.
- directionality_index and detect_tads: contact-matrix topology analysis.

## What the module does
The module unifies several related epigenomic workflows. It provides sparse coverage containers for signal tracks, peak summaries for ATAC/ChIP-style data, methylation experiments for CpG analysis, and 3D genome structures for contact-based analyses. Many of the functions transform raw sparse matrices or interval collections into derived results that can then be summarized or plotted.

## How the structs work together
SparseCoverageVector is the compact signal representation used by Epigenome. Peak and PeakSet represent the regions discovered from chromatin data. MethylationExperiment and MethylationResult handle locus-level methylation comparison. SingleCellChromatinExperiment wraps a single-cell RNA-like object with peak accessibility and reductions. ContactMatrix and TadResult support higher-order genome organization analysis.

## Typical usage
1. Create an Epigenome or PeakSet from interval and coverage data.
2. Use call_peaks and summarize_peak_support to identify and score peak regions.
3. Run differential_binding when comparing two conditions.
4. Build MethylationExperiment objects and run differential_methylation for locus-level comparisons.
5. For single-cell chromatin data, use tfidf, run_lsi, run_umap, cluster_cells, and calculate_coaccessibility.
6. Use directionality_index and detect_tads when working with contact matrices.

## Important implementation details
- The module reuses GenomicRanges interval logic so signals and regions stay coordinate-consistent.
- It also pulls in DifferentialExpression and SingleCell components so epigenomic analysis can share preprocessing and modeling code.
- DataFrame conversion methods are provided for peaks, support metrics, methylation results, and TAD results so output can be inspected easily.

## Why this file matters
This module is the package's epigenomics center. It keeps accessibility, methylation, single-cell chromatin, and contact-structure analysis in one coherent place instead of scattering those workflows across multiple unrelated modules.
