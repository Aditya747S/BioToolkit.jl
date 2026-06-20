# `search.jl` - Sequence Search, K-mer Indexing, and BLAST Parsing

## Overview

`search.jl` provides local sequence search utilities: compact k-mer indexing, BLAST-like seed-and-extend search, a wrapper for local BLAST+ commands, and lightweight parsers for BLAST XML and tabular output. It is designed to work with BioToolkit `BioSequence` objects while accepting string-like inputs for convenience.

### Purpose

Bioinformatics workflows often need two levels of search: quick in-process sequence lookup for small databases and interoperability with BLAST output from external tools. This module covers both without making HTTP/XML dependencies mandatory.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **K-mers are packed into `UInt64`** | Up to eight one-byte symbols fit into one integer key, making dictionary indexing simple and fast. |
| **Search accepts typed sequences first** | `BioSequence` dispatch preserves alphabet validation, while fallback methods coerce strings into inferred alphabets. |
| **Local search is seed-and-extend** | Exact k-mer matches identify candidate locations; ungapped X-drop and bounded gapped alignment refine them. |
| **Remote BLAST is intentionally a placeholder** | `qblast` requires optional HTTP/XML dependencies and raises a clear error instead of silently adding dependencies. |
| **Parsers are dependency-light** | BLAST XML parsing uses targeted string/regex extraction instead of requiring an XML package. |
| **File parsers are provenance-aware** | File reads can register hashes, source names, and record counts when provenance is active. |

---

## Table of Contents

1. [K-mer Index](#1-k-mer-index)
2. [Search Result Types](#2-search-result-types)
3. [Building an Index](#3-building-an-index)
4. [Local Seed-and-Extend Search](#4-local-seed-and-extend-search)
5. [Local BLAST+ Wrapper](#5-local-blast-wrapper)
6. [BLAST XML Types and Parsing](#6-blast-xml-types-and-parsing)
7. [BLAST Tabular Types and Parsing](#7-blast-tabular-types-and-parsing)
8. [Limitations](#8-limitations)
9. [Quick Reference](#9-quick-reference)

---

## 1. K-mer Index

### `KmerIndex`

```julia
struct KmerIndex
    k::Int
    database::Dict{UInt64, Vector{Tuple{Int, Int}}}
    target_names::Vector{String}
    target_sequences::Vector{Vector{UInt8}}
end
```

Stores a lookup table from packed k-mer to all target occurrences.

| Field | Description |
|---|---|
| `k` | K-mer width. Must be at most 8. |
| `database` | Maps packed k-mer keys to `(target_index, target_position)` hits. |
| `target_names` | Names parallel to input targets. |
| `target_sequences` | Copied raw byte data for each target. |

K-mers are packed by shifting one byte at a time into a `UInt64`.

---

## 2. Search Result Types

### `HighScoringPair`

```julia
struct HighScoringPair
    query_id::String
    target_id::String
    target_idx::Int
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    score::Float64
    evalue::Float64
    identity::Float64
    alignment_length::Int
end
```

Represents a local alignment hit. `run_blast` fills BLAST metadata from outfmt 6; `local_search` uses the legacy constructor for internally generated hits and leaves some BLAST-specific fields at default values.

---

## 3. Building an Index

### `build_index(targets, names=String[]; k=4)`

Builds a `KmerIndex` from target sequences.

Requirements:

- `k <= 8`
- Targets should be `BioSequence` values or string-like objects coercible to `BioSequence`.
- If `names` is empty, generated target names are used.

Example:

```julia
targets = AASeq.(["MTEYKLVVVG", "GAGGVGKSAL"])
index = build_index(targets, ["ras1", "ras2"]; k=3)
```

Complexity is O(total target length) for index construction.

---

## 4. Local Seed-and-Extend Search

### `local_search(query, index; scoring, x_drop=15, min_score=25, gap_open=-5, gap_extend=-1, envelope=20, use_threads=true, use_cuda=false)`

Performs a BLAST-like local search against a `KmerIndex`.

Algorithm:

1. Pack each query k-mer and look up exact seed matches.
2. For each seed, run ungapped X-drop extension.
3. Around promising seeds, run bounded gapped pairwise alignment in a local envelope.
4. Sort HSPs by descending score.
5. Filter overlapping query intervals greedily.

Important parameters:

| Parameter | Description |
|---|---|
| `scoring` | Required pairwise scoring object, usually a substitution matrix wrapper. |
| `x_drop` | Ungapped extension stops when score falls this far below the best extension score. |
| `min_score` | Minimum score required to keep an HSP. |
| `gap_open`, `gap_extend` | Gapped extension penalties. |
| `envelope` | Number of bases/residues around the seed extension used for bounded alignment. |
| `use_threads` | Uses Julia threads when multiple seeds are present. |
| `use_cuda` | Accepted by the API, but current logic still uses CPU/threaded seed processing. |

Example:

```julia
index = build_index(AASeq.(["MTEYKLVVVGAG", "GGGGKLVVVAAA"]); k=3)
scoring = MatrixPairwiseScoring(substitution_matrix("ACDEFGHIKLMNPQRSTVWY"))
hits = local_search(AASeq("KLVVVG"), index; scoring=scoring, min_score=10)
```

---

## 5. Local BLAST+ Wrapper

### `run_blast(program, query, database; options="", outfmt=6)`

Runs an installed BLAST+ executable through the system shell.

Behavior:

- Writes the query to a temporary FASTA file.
- Creates a temporary output file.
- Runs `program -query <query> -db <database> -out <output> -outfmt <outfmt> <options>`.
- Parses tabular outfmt 6 into `HighScoringPair` objects.
- Removes temporary files in a `finally` block.

Requirements:

- BLAST+ must be installed and available on `PATH`.
- The database must already be built for the chosen BLAST program.

### `qblast(program, database, query; options=Dict())`

Placeholder for remote BLAST. It throws an error explaining that optional `HTTP` and `EzXML` dependencies are required.

---

## 6. BLAST XML Types and Parsing

### `BlastXMLHSP{A,B}`

Stores one BLAST XML HSP with score/e-value statistics, coordinates, typed query and hit sequences, and midline text.

### `BlastXMLHit`

Stores a hit ID, description, hit length, and its HSPs.

### `BlastXMLRecord`

Stores one BLAST query result, including program/version/database metadata, query metadata, and hits.

### `parse_blast_xml(xml_or_io, provenance_source="parse_blast_xml", provenance_hash=nothing)`

Parses BLAST XML text or an `IO` object into `Vector{BlastXMLRecord}`.

The parser extracts:

- program, version, database
- query ID, definition, length
- hit ID, definition, length
- HSP bit score, e-value, identity count, positives, alignment length
- query and hit coordinates
- query sequence, hit sequence, midline

### `read_blast_xml(path)`

Reads a BLAST XML file and calls `parse_blast_xml`. With provenance enabled, file bytes are hashed before parsing.

---

## 7. BLAST Tabular Types and Parsing

### `BlastTabularHSP`

Stores coordinates, identity, alignment length, mismatch/gap counts, e-value, and bit score from BLAST outfmt 6.

### `BlastTabularHit`

Groups HSPs under a target ID.

### `BlastTabularRecord`

Groups hits under a query ID.

### `parse_blast_tabular(lines_or_io, provenance_source="parse_blast_tabular", provenance_hash=nothing)`

Parses BLAST outfmt 6 tabular text. Comment and empty lines are ignored. Lines with fewer than twelve fields are skipped.

Expected field order:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

### `read_blast_tabular(path)`

Reads a BLAST tabular file and parses it into query-grouped records.

---

## 8. Limitations

- `qblast` is not implemented in this dependency-light core module.
- XML parsing is targeted to standard BLAST XML tags and is not a general XML parser.
- `local_search` uses exact k-mer seeds; it will miss alignments without an exact seed of width `k`.
- `build_index` uses one byte per symbol and requires `k <= 8`.

---

## 9. Quick Reference

| Function/Type | Purpose |
|---|---|
| `KmerIndex` | Packed k-mer lookup table. |
| `HighScoringPair` | Local search or BLAST hit. |
| `build_index` | Build a k-mer database from targets. |
| `local_search` | BLAST-like seed-and-extend search. |
| `run_blast` | Run local BLAST+ and parse outfmt 6. |
| `qblast` | Remote BLAST placeholder. |
| `BlastXMLRecord` | Parsed BLAST XML query result. |
| `BlastXMLHit` | Parsed BLAST XML hit. |
| `BlastXMLHSP` | Parsed BLAST XML HSP. |
| `parse_blast_xml` | Parse XML text or IO. |
| `read_blast_xml` | Read and parse XML file. |
| `BlastTabularRecord` | Parsed tabular query result. |
| `BlastTabularHit` | Parsed tabular hit. |
| `BlastTabularHSP` | Parsed tabular HSP. |
| `parse_blast_tabular` | Parse outfmt 6 text or IO. |
| `read_blast_tabular` | Read and parse outfmt 6 file. |

