# `genomicranges.jl` - Genomic Interval Operations

## Overview

`genomicranges.jl` defines the `BioToolkit.GenomicRanges` submodule. It provides typed genomic intervals, interval collections, efficient overlap queries, nearest-neighbor operations, interval transforms, set operations, coverage generation, and table conversion.

### Purpose

Genomic analysis relies heavily on interval algebra: overlap peaks with genes, trim ranges to chromosome lengths, calculate coverage, derive promoters, disjoin intervals, and subtract masked regions. This module provides those operations in pure Julia using sorted interval collections with per-chromosome indices and interval trees.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Closed interval coordinates** | `GenomicInterval(left, right)` stores inclusive endpoints. |
| **Chromosome partitioning** | `IntervalCollection` stores per-chromosome ranges for faster queries. |
| **Sorted arrays plus interval trees** | Collections keep sorted arrays and per-chromosome `IntervalTree{Int}` indices. |
| **Metadata carries provenance** | Interval metadata is stamped unless provenance already exists. |
| **Set operations preserve interval metadata** | Generated intervals usually copy metadata from the source interval. |
| **DataFrame conversion preserves metadata columns** | Interval metadata keys become extra columns in tabular output. |

---

## Table of Contents

1. [Core Types](#1-core-types)
2. [Constructors](#2-constructors)
3. [Collection Building](#3-collection-building)
4. [Overlap and Distance Queries](#4-overlap-and-distance-queries)
5. [Interval Transforms](#5-interval-transforms)
6. [Set Operations](#6-set-operations)
7. [Genome Bounds and Gaps](#7-genome-bounds-and-gaps)
8. [Coverage](#8-coverage)
9. [DataFrame and File I/O](#9-dataframe-and-file-io)
10. [Quick Reference](#10-quick-reference)

---

## 1. Core Types

### `GenomicInterval`

```julia
struct GenomicInterval
    chrom::String
    left::Int
    right::Int
    strand::Char
    metadata::Dict{String,Any}
end
```

**Kind:** Concrete struct

**Description:** One closed genomic interval on one chromosome.

**Validation:**

- strand must be `'+'`, `'-'`, or `'.'`;
- endpoints are normalized with `normalize_interval`;
- metadata keys are converted to strings;
- metadata is provenance-stamped when missing.

### `IntervalCollection`

```julia
struct IntervalCollection
    chroms::PooledStringVector
    starts::Vector{Int}
    ends::Vector{Int}
    strands::Vector{Char}
    metadata::Vector{Dict{String,Any}}
    intervals::Vector{GenomicInterval}
    chrom_indices::Dict{String,UnitRange{Int}}
    chrom_end_indices::Dict{String,Vector{Int}}
    trees::Dict{String,IntervalTree{Int}}
end
```

**Description:** Sorted, indexed collection of genomic intervals. It stores both columnar arrays and the original interval vector.

### `CoverageSegment`

```julia
struct CoverageSegment
    chrom::String
    start::Int
    stop::Int
    depth::Int
end
```

**Description:** Constant-depth coverage segment over a closed coordinate range.

### `SeqInfo`

```julia
struct SeqInfo
    chrom::String
    length::Int
    is_circular::Bool
end
```

**Description:** Chromosome metadata used by trim, gaps, and complement operations.

---

## 2. Constructors

### `GenomicInterval`

```julia
GenomicInterval(chrom, left, right)
GenomicInterval(chrom, left, right, strand)
GenomicInterval(chrom, left, right, strand, metadata)
GenomicInterval(record::BedRecord)
GenomicInterval(record::GffRecord)
```

**Description:** Creates intervals from explicit coordinates or parsed BED/GFF records. BED records convert 0-based BED starts to 1-based closed starts with `start + 1`.

**Example:**

```julia
interval = GenomicInterval("chr1", 100, 200, '+', Dict("gene" => "TP53"))
```

---

## 3. Collection Building

### `build_collection`

```julia
build_collection(intervals::AbstractVector{<:GenomicInterval}) -> IntervalCollection
IntervalCollection(intervals::AbstractVector{<:GenomicInterval})
```

**Description:** Sorts intervals by chromosome, start, end, and strand; builds pooled chromosome arrays; creates per-chromosome index ranges; and builds interval trees for overlap queries.

**Interface:**

```julia
length(collection)
isempty(collection)
iterate(collection)
collection[i]
```

---

## 4. Overlap and Distance Queries

### `find_overlaps` / `overlap`

```julia
find_overlaps(query::GenomicInterval, subject::IntervalCollection)
overlap(query, subject)
```

**Description:** Returns subject intervals overlapping the query. The collection uses chromosome-specific interval trees for candidate lookup.

### `distance`

```julia
distance(query::GenomicInterval, subject::IntervalCollection)
```

**Description:** Computes distances from a query to intervals in the subject collection.

### `nearest` / `find_nearest`

```julia
nearest(query, subject; select=:all)
find_nearest(query, subject; select=:all)
```

**Description:** Finds nearest intervals. `select` controls whether all nearest ties or a subset are returned.

### `precede` and `follow`

```julia
precede(query, subject)
follow(query, subject)
```

**Description:** Finds the nearest downstream or upstream interval relative to the query.

### `find_overlaps_parallel`

```julia
find_overlaps_parallel(query::IntervalCollection, subject::IntervalCollection; multi_thread=true)
```

**Description:** Computes overlaps for a collection of query intervals, optionally using Julia threads.

---

## 5. Interval Transforms

### `shift`

```julia
shift(intervals, delta)
shift(collection, delta)
```

**Description:** Adds `delta` to interval starts and ends.

### `flank`

```julia
flank(intervals, width; start=true)
```

**Description:** Creates upstream or downstream flanking intervals, respecting strand orientation.

### `resize`

```julia
resize(intervals, width; fix=:start)
```

**Description:** Resizes intervals to `width`, anchoring at `:start`, `:center`, or `:end`.

### `promoters`

```julia
promoters(intervals, upstream, downstream)
```

**Description:** Constructs promoter-like intervals around transcription starts, respecting strand.

### `narrow`

```julia
narrow(intervals, start=0, stop=0)
```

**Description:** Trims bases from interval starts and ends. Intervals that invert are skipped.

---

## 6. Set Operations

### `setdiff`

```julia
setdiff(left::GenomicInterval, right::IntervalCollection)
setdiff(left::IntervalCollection, right::IntervalCollection)
```

**Description:** Subtracts blocking intervals using sweep-line interval subtraction.

### `reduce`

```julia
reduce(collection::IntervalCollection)
```

**Description:** Merges overlapping intervals within the collection.

### `intersect` and `union`

```julia
intersect(left::IntervalCollection, right::IntervalCollection)
union(left::IntervalCollection, right::IntervalCollection)
```

**Description:** Computes interval intersections or merged union collections.

### Pairwise operations

```julia
pintersect(x, y)
punion(x, y)
psetdiff(x, y)
```

**Description:** Pairwise interval operations over equally sized interval vectors.

---

## 7. Genome Bounds and Gaps

### `trim`

```julia
trim(intervals, seqlengths::Dict{String,Int})
trim(collection, seqlengths)
trim(intervals, seqinfo::Dict{String,SeqInfo})
```

**Description:** Clips intervals to chromosome bounds.

### `gaps`

```julia
gaps(intervals, seqlengths)
gaps(collection, seqlengths)
```

**Description:** Returns uncovered intervals for each chromosome.

### `complement`

```julia
complement(intervals, seqlengths)
complement(collection, seqlengths)
```

**Description:** Alias-style wrapper around `gaps`.

### `disjoin`

```julia
disjoin(collection::IntervalCollection)
```

**Description:** Splits overlapping intervals into disjoint atomic pieces.

---

## 8. Coverage

### `coverage`

```julia
coverage(collection::IntervalCollection) -> Vector{CoverageSegment}
coverage(intervals::AbstractVector{<:GenomicInterval}) -> Vector{CoverageSegment}
```

**Description:** Computes constant-depth coverage segments across intervals.

---

## 9. DataFrame and File I/O

### `DataFrame(intervals)`

```julia
DataFrames.DataFrame(intervals::AbstractVector{<:GenomicInterval})
```

**Description:** Converts intervals into a DataFrame with core coordinate columns plus metadata-derived columns.

### `read_intervals`

```julia
read_intervals(df::DataFrame)
read_intervals(io::IO)
read_intervals(path::String)
```

**Description:** Reads intervals from DataFrame or delimited text. Text input expects a header row and coordinate columns.

---

## 10. Quick Reference

| API | Purpose |
|---|---|
| `GenomicInterval` | One closed genomic interval. |
| `IntervalCollection` | Sorted indexed interval set. |
| `CoverageSegment` | Constant-depth coverage run. |
| `SeqInfo` | Chromosome length/circularity metadata. |
| `build_collection` | Build indexed collection. |
| `find_overlaps` | Query overlaps. |
| `nearest` | Find nearest intervals. |
| `precede` / `follow` | Downstream/upstream neighbors. |
| `shift`, `flank`, `resize` | Coordinate transforms. |
| `promoters`, `narrow` | Common genomic range transforms. |
| `setdiff`, `intersect`, `union` | Set operations. |
| `trim`, `gaps`, `complement` | Genome-bound operations. |
| `disjoin` | Split into non-overlapping pieces. |
| `coverage` | Coverage segments. |
| `read_intervals` | Read interval tables. |

---

## Complete Usage Example

```julia
using BioToolkit

intervals = [
    GenomicInterval("chr1", 100, 200, '+'),
    GenomicInterval("chr1", 150, 250, '+'),
    GenomicInterval("chr2", 10, 50, '.'),
]

collection = build_collection(intervals)
query = GenomicInterval("chr1", 180, 220)

hits = find_overlaps(query, collection)
merged = reduce(collection)
prom = promoters(intervals, 1000, 200)
cov = coverage(collection)
```
