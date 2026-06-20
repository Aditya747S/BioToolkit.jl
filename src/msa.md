# `msa.jl` - Multiple Sequence Alignment I/O and Manipulation

## Overview

`msa.jl` provides BioToolkit's native multiple sequence alignment container plus readers, writers, slicing, consensus lines, column statistics, and alignment concatenation. It is designed around `SeqRecordLite` and typed `BioSequence` values, so parsed alignments remain compatible with the rest of the toolkit's sequence, motif, phylogeny, and provenance APIs.

### Purpose

Multiple sequence alignments appear in many incompatible text formats. This module gives BioToolkit one common in-memory representation and enough format support to move between FASTA, CLUSTAL, Stockholm, PIR, NEXUS, MSF, PHYLIP, MAF, and GCG-style data without immediately depending on external alignment packages.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Alignment rows are `SeqRecordLite` values** | Row identifiers, descriptions, annotations, and letter annotations are preserved instead of storing raw sequence strings only. |
| **All rows must share one alignment width** | The constructor validates row lengths so downstream column operations never have to handle ragged alignments. |
| **Column annotations are first-class** | Formats such as Stockholm and CLUSTAL carry consensus or per-column annotation lines; those values are stored in `column_annotations`. |
| **Format detection is conservative** | `read_alignment(..., "auto")` identifies obvious signatures and falls back to FASTA when no stronger signal is present. |
| **Slicing preserves metadata where possible** | Row and column slicing returns new alignment objects and slices compatible column annotations. |
| **Provenance-aware file reads** | When a provenance context is active, file bytes are hashed and each alignment plus each parsed record is registered. |

---

## Table of Contents

1. [Core Types](#1-core-types)
2. [Constructors](#2-constructors)
3. [Reading Alignments](#3-reading-alignments)
4. [Container Interface](#4-container-interface)
5. [Slicing and Mutation](#5-slicing-and-mutation)
6. [Consensus and Column Statistics](#6-consensus-and-column-statistics)
7. [Writing Alignments](#7-writing-alignments)
8. [Format-Specific Notes](#8-format-specific-notes)
9. [Quick Reference](#9-quick-reference)

---

## 1. Core Types

### `AbstractMultipleSequenceAlignment`

```julia
abstract type AbstractMultipleSequenceAlignment end
```

Abstract supertype for alignment containers. `MultipleSequenceAlignment` is the concrete implementation currently provided by BioToolkit.

### `MultipleSequenceAlignment`

```julia
mutable struct MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    records::Vector{SeqRecordLite}
    annotations::Dict{Symbol,Any}
    column_annotations::Dict{Symbol,Any}
end
```

Stores a row-oriented multiple sequence alignment.

| Field | Type | Description |
|---|---|---|
| `records` | `Vector{SeqRecordLite}` | One aligned sequence per row. Every sequence must have the same length. |
| `annotations` | `Dict{Symbol,Any}` | Alignment-level metadata. |
| `column_annotations` | `Dict{Symbol,Any}` | Per-column annotations such as CLUSTAL consensus or Stockholm `#=GC` tracks. |

`column_annotations` values may be strings or vectors. If a value has a length, that length must match the alignment width.

---

## 2. Constructors

### `MultipleSequenceAlignment(records; annotations, column_annotations)`

```julia
MultipleSequenceAlignment(
    records::AbstractVector;
    annotations::AbstractDict=Dict{Symbol,Any}(),
    column_annotations::AbstractDict=Dict{Symbol,Any}())
```

Builds an alignment from strings or `SeqRecordLite` objects.

Behavior:

- `SeqRecordLite` rows are copied into new records.
- `String` rows are converted to typed `BioSequence` values using inferred alphabet.
- Empty input is allowed and produces a zero-row, zero-column alignment.
- Nonempty input must have identical sequence lengths.
- Column annotations with measurable length must match the alignment width.

Example:

```julia
aln = MultipleSequenceAlignment(["ACGT", "A-GT", "AC-T"])
size(aln)                 # (3, 4)
get_alignment_length(aln) # 4
```

---

## 3. Reading Alignments

### `read_alignment(source::String, format="auto")`

Reads an alignment from a file path.

```julia
alignment = read_alignment("family.sto", "stockholm")
```

When provenance is active, the file is read as bytes first, SHA-256 hashed, and then parsed from an `IOBuffer`.

### `read_alignment(io::IO, format="auto")`

Reads an alignment from any Julia `IO` object.

Supported formats:

| Format argument | Reader |
|---|---|
| `"auto"` | Detects format from the first line. |
| `"fasta"`, `"fa"` | FASTA alignment. |
| `"clustal"` | CLUSTAL/MUSCLE block alignment. |
| `"stockholm"`, `"sto"` | Stockholm 1.0. |
| `"pir"` | PIR/NBRF-style alignment. |
| `"nexus"` | NEXUS data matrix. |
| `"msf"` | GCG/MSF alignment. |
| `"phylip"`, `"phylip-relaxed"` | Sequential or interleaved PHYLIP-like alignment. |
| `"maf"` | Multiple Alignment Format blocks. |
| `"gcg"` | GCG alignment flavor. |

Unsupported formats throw `ArgumentError`.

---

## 4. Container Interface

`MultipleSequenceAlignment` behaves like a row container with matrix-like indexing.

| Method | Description |
|---|---|
| `length(alignment)` | Number of rows. |
| `get_alignment_length(alignment)` | Number of columns. |
| `size(alignment)` | `(rows, columns)`. |
| `iterate(alignment)` | Iterates over `SeqRecordLite` rows. |
| `copy(alignment)` | Returns a new alignment object with copied row records. |
| `alignment[i]` | Returns row `i` as a `SeqRecordLite`. |
| `alignment[:, j]` | Returns column `j` as a string. |
| `alignment[i, j]` | Returns a single aligned symbol. |
| `alignment[i, :]` | Returns row `i`. |
| `alignment[:, cols]` | Returns a column-sliced alignment. |
| `alignment[rows, :]` | Returns a row-sliced alignment. |
| `alignment[rows, cols]` | Returns a row-and-column sliced alignment. |

Example:

```julia
aln = MultipleSequenceAlignment(["ACGT", "A-GT"])
aln[:, 2]   # "C-"
aln[1, 3]   # 'G'
aln[:, 2:4] # alignment with 3 columns
```

---

## 5. Slicing and Mutation

### `push!(alignment, record)`

Adds a `SeqRecordLite` row. If the alignment is nonempty, the new row length must equal the existing alignment width.

### `append!(alignment, record)`

Alias for `push!`.

### `deleteat!(alignment, index)`

Deletes one or more rows from the alignment.

### `sort!(alignment; key=nothing, reverse=false)`

Sorts rows in place. By default rows are sorted by `record.identifier`. A custom `key` function can be provided.

### `left + right`

Concatenates two alignments column-wise.

Requirements and behavior:

- Both alignments must have the same number of rows.
- Rows are paired by position, not by identifier.
- Only identical alignment-level annotation keys are preserved.
- Shared column annotation keys are concatenated when possible.

---

## 6. Consensus and Column Statistics

### `alignment_column_counts(alignment, col; gap='-')`

Returns `Dict{Char,Int}` counts for non-gap symbols in one column.

### `alignment_column_frequencies(alignment, col; gap='-')`

Returns non-gap symbol frequencies for one column. Fully gapped columns return an empty dictionary.

### `alignment_column_symbol(alignment, col; scoring=nothing, strong_threshold=1, weak_threshold=0, style=:clustal)`

Computes the CLUSTAL/EMBOSS-style consensus character for a single column.

Symbols:

| Output | Meaning |
|---|---|
| `*` | All non-gap symbols are identical in CLUSTAL style. |
| `|` | All non-gap symbols are identical in EMBOSS style. |
| `:` | All pairwise scores meet the strong threshold. |
| `.` | At least one pairwise score meets the weak threshold. |
| space | No consensus class. |

### `alignment_symbol_line(alignment; kwargs...)`

Computes the consensus/symbol line across all columns. If `column_annotations[:clustal_consensus]` exists and `style == :clustal`, that stored annotation is returned directly.

### `clustal_consensus(alignment; kwargs...)`

Convenience wrapper for `alignment_symbol_line(..., style=:clustal)`.

### `emboss_consensus(alignment; kwargs...)`

Convenience wrapper for `alignment_symbol_line(..., style=:emboss)`.

---

## 7. Writing Alignments

The module contains format writers used by the public alignment output API.

Supported output formats include:

| Format | Notes |
|---|---|
| FASTA | Writes one row per FASTA record, with configurable wrapping. |
| Stockholm | Writes `# STOCKHOLM 1.0` plus compatible column annotations. |
| MSF | Writes GCG/MSF-style header and 10-character grouped blocks. |
| PIR | Writes `>P1;identifier` records ending in `*`. |
| NEXUS | Writes a simple `Begin data` matrix. |
| PHYLIP | Writes strict 10-character identifiers unless relaxed mode is requested. |

Empty alignments cannot be formatted and throw `ArgumentError`.

---

## 8. Format-Specific Notes

### CLUSTAL

CLUSTAL and MUSCLE headers are skipped. Consensus lines containing only `*`, `:`, `.`, `|`, and spaces are collected into `column_annotations[:clustal_consensus]`.

### Stockholm

`#=GC` lines are parsed into `column_annotations`. Other comments are ignored.

### PHYLIP

The reader accepts a header of sequence count and alignment length, then supports continuation blocks. Parsed record counts and sequence lengths are validated.

### PIR

The parser expects headers containing a semicolon, skips the description line, strips whitespace from sequence lines, and removes terminal `*` markers.

### NEXUS

The parser reads `DIMENSIONS` for `ntax` and `nchar` when present and validates both values against parsed rows.

---

## 9. Quick Reference

| Function | Returns | Notes |
|---|---|---|
| `read_alignment(path, format)` | `MultipleSequenceAlignment` | File parser with optional provenance. |
| `read_alignment(io, format)` | `MultipleSequenceAlignment` | IO parser. |
| `MultipleSequenceAlignment(records)` | `MultipleSequenceAlignment` | Validates equal row lengths. |
| `get_alignment_length(aln)` | `Int` | Alignment width. |
| `alignment_column_counts(aln, col)` | `Dict{Char,Int}` | Excludes gaps. |
| `alignment_column_frequencies(aln, col)` | `Dict{Char,Float64}` | Excludes gaps. |
| `alignment_symbol_line(aln)` | `String` | CLUSTAL/EMBOSS consensus logic. |
| `clustal_consensus(aln)` | `String` | CLUSTAL wrapper. |
| `emboss_consensus(aln)` | `String` | EMBOSS wrapper. |
| `push!(aln, record)` | `MultipleSequenceAlignment` | Adds a row. |
| `sort!(aln)` | `MultipleSequenceAlignment` | Sorts rows in place. |
| `left + right` | `MultipleSequenceAlignment` | Concatenates columns. |

