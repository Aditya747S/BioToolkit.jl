# `query.jl` - Genomic Table Queries and Coverage Utilities

## Overview

`query.jl` provides table-level genomic query helpers: region filtering, position binning, interval-window coverage, BigWig writing, and basic interval arithmetic.

### Purpose

Genomic workflows repeatedly need the same operations:

- filter a table to one chromosome/range;
- count positions in bins;
- aggregate interval coverage by fixed windows;
- write coverage tracks;
- normalize, intersect, merge, and subtract intervals.

This file implements those operations over generic Tables.jl-compatible inputs where possible and uses optional CUDA acceleration for selected hot paths.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Tables.jl column access** | `filter_region` works with any table that supports `Tables.columntable`. |
| **Sorted and unsorted paths** | Sorted region queries use binary search; unsorted queries use masks. |
| **Threaded CPU implementations** | Position binning and window coverage can use per-thread partial dictionaries. |
| **Optional CUDA path** | GPU kernels are loaded lazily through `lazy_gpu.jl`. |
| **Half-open interval logic for windows** | Window coverage treats intervals as `[start, stop)` during overlap accumulation. |
| **Provenance-aware results** | Public query helpers stamp or register provenance when a context is active. |

---

## 1. Region Filtering

### `filter_region`

```julia
filter_region(table, chrom, min_pos, max_pos; sorted=false, prov_ctx=nothing)
```

**Description:** Filters rows in a genomic table by chromosome and coordinate range.

**Expected columns:** `chrom`, `pos`, `id`, `ref`, `alt`, `qual`.

**Parameters:**

| Parameter | Description |
|---|---|
| `table` | Tables.jl-compatible table. |
| `chrom` | Chromosome string to match. |
| `min_pos` | Minimum position, inclusive. |
| `max_pos` | Maximum position, inclusive. |
| `sorted` | Use sorted binary-search path when table is sorted by chromosome/position. |

**Returns:** Named tuple of filtered columns.

**Implementation paths:**

| Helper | Use |
|---|---|
| `_filter_region_mask` | Boolean mask over unsorted columns. |
| `_filter_region_sorted` | Binary search over sorted chromosome/position columns. |

---

## 2. Position Binning

### `bin_positions`

```julia
bin_positions(positions, bin_size; threaded=true, use_cuda=false, prov_ctx=nothing)
```

**Description:** Counts integer positions into zero-based bins using `fld(position, bin_size)`.

**Backends:**

| Backend | Helper |
|---|---|
| Serial CPU | `_bin_positions_serial` |
| Threaded CPU | `_bin_positions_threaded` |
| CUDA | `_CUDA_QUERY_BIN_IMPL` after `_ensure_cuda_query!()` |

**Returns:** `Dict{Int,Int}` mapping bin index to count.

**Example:**

```julia
bin_positions([10, 11, 25, 40], 10)
# Dict(1=>2, 2=>1, 4=>1)
```

---

## 3. Coverage Histograms

### `coverage_histogram`

```julia
coverage_histogram(table, chrom, bin_size; threaded=true, use_cuda=false, prov_ctx=nothing)
```

**Description:** Filters a genomic table to `chrom` and bins its `pos` column.

**Returns:** `Dict{Int,Int}` of bin counts.

**Implementation:** Calls `filter_region(table, chrom, typemin(Int), typemax(Int))`, then `bin_positions`.

---

## 4. Window Coverage

### `window_coverage(starts, stops, window_size)`

```julia
window_coverage(starts, stops, window_size; threaded=true, use_cuda=false, prov_ctx=nothing)
```

**Description:** Computes total interval overlap per fixed-width window.

**Validation:**

- `window_size` must be positive.
- `length(starts) == length(stops)`.

**Returns:** `Dict{Int,Int}` mapping zero-based window index to covered bases.

### `window_coverage(table, chrom, window_size)`

```julia
window_coverage(table, chrom, window_size; sorted=false, threaded=true, use_cuda=false)
```

**Description:** Filters a genomic interval table by chromosome and computes window coverage from start/stop columns.

---

## 5. BigWig Writing

### `write_bigwig(coverage_vector, output_path)`

```julia
write_bigwig(coverage_vector::AbstractVector{<:Real}, output_path::String; chrom="chr1", start=1)
```

**Description:** Writes a single-chromosome coverage vector to BigWig. Consecutive equal values are collapsed into runs.

### `write_bigwig(coverage_vectors, output_path)`

```julia
write_bigwig(coverage_vectors::AbstractDict, output_path::String; starts=1)
```

**Description:** Writes a dictionary of chromosome => coverage vector to BigWig. `starts` can be a scalar or a per-chromosome mapping.

**Internal helpers:**

| Helper | Purpose |
|---|---|
| `_write_bigwig_runs` | Convert coverage vector into run records. |
| `_bigwig_start` | Resolve per-chromosome start position. |
| `_write_bigwig_coverage_map` | Write multi-chromosome coverage map. |

---

## 6. Interval Arithmetic

### `normalize_interval`

```julia
normalize_interval(start, stop) -> Tuple{Int,Int}
normalize_interval(interval::Tuple) -> Tuple{Int,Int}
```

**Description:** Returns an ordered `(start, stop)` tuple regardless of input order.

### `interval_length`

```julia
interval_length(start, stop) -> Int
```

**Description:** Length of a closed interval after normalization.

### `interval_contains`

```julia
interval_contains(start, stop, position) -> Bool
```

**Description:** Tests whether a position lies inside a closed interval.

### `interval_overlaps`

```julia
interval_overlaps(left_start, left_stop, right_start, right_stop) -> Bool
```

**Description:** Tests whether two closed intervals overlap.

### `interval_intersection`

```julia
interval_intersection(left, right) -> Union{Tuple{Int,Int},Nothing}
```

**Description:** Returns the closed interval intersection or `nothing`.

### `interval_union`

```julia
interval_union(left, right) -> Tuple{Int,Int}
```

**Description:** Returns the minimum spanning interval covering both intervals.

### `merge_intervals`

```julia
merge_intervals(intervals) -> Vector{Tuple{Int,Int}}
```

**Description:** Sorts and merges overlapping or adjacent closed intervals.

### `interval_difference`

```julia
interval_difference(left, right) -> Vector{Tuple{Int,Int}}
```

**Description:** Subtracts one interval from another and returns remaining pieces.

---

## Quick Reference

| API | Purpose |
|---|---|
| `filter_region` | Filter genomic table by chromosome/range. |
| `bin_positions` | Count positions into bins. |
| `coverage_histogram` | Build position histogram for a chromosome. |
| `window_coverage` | Aggregate interval coverage by windows. |
| `write_bigwig` | Write coverage tracks. |
| `normalize_interval` | Sort interval endpoints. |
| `interval_length` | Closed interval length. |
| `interval_contains` | Position membership. |
| `interval_overlaps` | Interval overlap test. |
| `interval_intersection` | Shared interval region. |
| `interval_union` | Spanning interval. |
| `merge_intervals` | Merge interval list. |
| `interval_difference` | Subtract interval. |

---

## Complete Usage Example

```julia
table = (
    chrom = ["chr1", "chr1", "chr2"],
    pos = [100, 250, 100],
    id = ["a", "b", "c"],
    ref = ["A", "G", "T"],
    alt = ["C", "A", "G"],
    qual = [60.0, 50.0, 40.0],
)

region = filter_region(table, "chr1", 1, 200)
bins = bin_positions(region.pos, 100)

starts = [0, 50, 120]
stops = [100, 180, 200]
cov = window_coverage(starts, stops, 50)

interval_overlaps(10, 20, 18, 25)
merge_intervals([(1, 5), (4, 10), (20, 25)])
```
