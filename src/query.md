# query.jl

## Purpose
This file implements fast region filtering and binning for genomic coordinate data. It is an optimization-oriented utility layer used by higher-level interval and coverage workflows.

## Main functions
- `filter_region(table, chrom, min_pos, max_pos; sorted=false)` filters table-like input to a genomic interval.
- `bin_positions(positions, bin_size; threaded=false)` bins positions into histogram buckets.
- `coverage_histogram(table, chrom, bin_size; threaded=false)` filters a chromosome and then bins positions.
- `window_coverage(starts, stops, window_size)` computes coverage over fixed-width windows.
- `window_coverage(table, chrom, window_size; sorted=false)` performs the same operation on table input.
- `write_bigwig(coverage_vector, output_path; chrom, start)` writes a BigWig track from one coverage vector.
- `write_bigwig(coverage_vectors, output_path; starts)` writes a multi-chromosome BigWig file.
- `normalize_interval`, `interval_length`, `interval_contains`, `interval_overlaps`, `interval_intersection`, `interval_union`, `merge_intervals`, and `interval_difference` provide basic interval arithmetic.

## How it is used
The normal usage pattern is to narrow a large coordinate table to a chromosome or range, then summarize it with `bin_positions` or `window_coverage`. The interval helpers support the rest of the package when regions need to be compared, merged, or clipped.

`write_bigwig` is the output path when users want to turn coverage vectors into browser-friendly tracks.

## Implementation notes
- The module supports both masked and sorted-table filtering strategies.
- GPU-aware query paths are routed through the lazy CUDA helper module when `use_cuda=true` is passed.
- The binning and coverage helpers promote CPU vectors to CUDA arrays on demand, so the opt-in GPU path works from normal table or array inputs.
- The interval helpers use half-open style comparisons, which keeps overlap logic consistent.

## Why it matters
Large genomic tables are only useful if they can be queried efficiently. This file keeps the filtering, binning, and track-writing logic centralized so higher-level modules can stay focused on biology rather than bookkeeping.
