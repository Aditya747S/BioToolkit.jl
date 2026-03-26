# genomicranges.jl

## Purpose
This file is the genomic interval engine for BioToolkit. It defines the coordinate model used by the package and provides the query and transformation operations needed to work with genomic annotations, coverage, and feature arithmetic.

## Main structs
- GenomicInterval: a chromosome interval with left and right coordinates, strand, and arbitrary metadata.
- IntervalCollection: a sorted, indexed set of genomic intervals with chromosome lookup tables.
- CoverageSegment: a depth segment used for binned coverage summaries.
- SeqInfo: metadata for a sequence or chromosome, including length and circularity.

## Public constructors and functions
- GenomicInterval(chrom, left, right): create an unstranded interval.
- GenomicInterval(chrom, left, right, strand): create a stranded interval.
- GenomicInterval(record::BedRecord): convert a BED record into a 1-based interval.
- GenomicInterval(record::GffRecord): convert a GFF record into a genomic interval.
- build_collection(intervals): sort intervals and build an indexed IntervalCollection.
- IntervalCollection(intervals): convenience constructor that delegates to build_collection.
- read_intervals: load interval-like records into genomic intervals.
- overlap and find_overlaps: test and search for overlapping features.
- nearest and find_nearest: locate the closest feature to a coordinate or interval.
- follow and precede: find adjacent intervals in genomic order.
- shift, flank, resize, promoters, narrow: generate derived intervals.
- trim, gaps, complement, disjoin, pintersect, punion, psetdiff: interval set and coverage-style operations.
- coverage: compute coverage summaries from interval collections.

## What the module does
The file turns raw genomic coordinates into a structured and searchable representation. Once intervals are collected, they are sorted and indexed by chromosome, which lets the package answer overlap and nearest-neighbor questions efficiently. The same data structure also supports common interval transformations such as shifting regions, resizing windows, generating promoter spans, and narrowing intervals.

## How the structs work together
GenomicInterval is the primary building block. IntervalCollection is the faster, indexed form that stores the same intervals in sorted order and provides per-chromosome access tables. CoverageSegment represents a summarized signal interval. SeqInfo stores sequence metadata that can be used when a caller needs context for coordinate-based operations.

## Typical usage
1. Create GenomicInterval objects manually or by converting BED/GFF records.
2. Call build_collection to obtain a searchable indexed collection.
3. Use overlap or find_overlaps to identify matching features.
4. Use nearest or find_nearest when you need the closest annotated region.
5. Apply shift, resize, promoters, or narrow to build derived intervals for plotting or downstream analysis.

## Important implementation details
- Intervals are normalized with normalize_interval so coordinates are stored consistently.
- build_collection sorts by chromosome and coordinate, then creates chromosome index ranges and end-index ranges.
- The collection implements length, isempty, iterate, and indexing so it behaves like a normal Julia container.
- Private helpers such as _find_last_start_leq support fast binary-search style querying.

## Why this file matters
This module is the coordinate backbone of BioToolkit. Nearly every genome-aware analysis depends on it, so keeping this logic centralized prevents duplicated interval math across the package.
