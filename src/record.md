# `record.jl` - Biological Sequence Record Types

## Overview

`record.jl` defines typed biological sequence records that bundle `BioSequence` objects with identifiers, descriptions, metadata, and quality scores.

### Purpose

Raw `BioSequence` values are ideal for algorithms, but real FASTA, FASTQ, and annotated sequence workflows need identifiers, descriptions, annotations, and per-letter metadata. `record.jl` provides those record containers while preserving BioToolkit's alphabet type safety.

Records are parametric on `BioAlphabet`, so a `FastqRecord{DNAAlphabet}` remains distinct from a protein record at the type level.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Records are alphabet-parametric** | `SeqRecord{A}`, `FastqRecord{A}`, and `SeqRecordLite{A}` carry the alphabet in the type. |
| **Metadata carries provenance** | Constructors ensure a provenance id in metadata or annotations dictionaries. |
| **FASTQ quality is stored as raw ASCII** | FASTQ quality strings are retained in Phred+33 form; integer conversion is explicit via `quality_scores`. |
| **Lightweight compatibility record** | `SeqRecordLite` supports MSA, annotation, and quality modules with Biopython-like fields. |
| **Equality ignores provenance metadata** | Two biologically identical records compare equal even if their provenance ids differ. |

---

## Table of Contents

1. [SeqRecord](#1-seqrecord)
2. [FastqRecord](#2-fastqrecord)
3. [SeqRecordLite](#3-seqrecordlite)
4. [Constructors](#4-constructors)
5. [Utility Methods](#5-utility-methods)
6. [Quality Scores](#6-quality-scores)
7. [Quick Reference](#7-quick-reference)

---

## 1. `SeqRecord`

```julia
struct SeqRecord{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    metadata::Dict{Symbol, Any}
end
```

**Kind:** Parametric concrete struct

**Description:** Annotated sequence record for FASTA-like and general sequence workflows.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `sequence` | `BioSequence{A}` | Typed biological sequence. |
| `identifier` | `String` | Record id. |
| `description` | `String` | Human-readable description. |
| `metadata` | `Dict{Symbol,Any}` | Record metadata and provenance. |

---

## 2. `FastqRecord`

```julia
struct FastqRecord{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    quality::String
    metadata::Dict{Symbol, Any}
end
```

**Kind:** Parametric concrete struct

**Description:** FASTQ record pairing a typed biological sequence with a raw Phred+33 ASCII quality string.

**Invariant:** `length(sequence) == length(quality)`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `sequence` | `BioSequence{A}` | Typed sequence. |
| `identifier` | `String` | FASTQ id. |
| `description` | `String` | FASTQ description. |
| `quality` | `String` | Raw Phred+33 quality string. |
| `metadata` | `Dict{Symbol,Any}` | Record metadata and provenance. |

---

## 3. `SeqRecordLite`

```julia
struct SeqRecordLite{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    name::String
    description::String
    annotations::Dict{Symbol, Any}
    letter_annotations::Dict{Symbol, Any}
end
```

**Kind:** Parametric concrete struct

**Description:** Lightweight record compatible with modules that need both record-level annotations and per-letter annotations. It is useful for MSA, annotation, and quality workflows.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `sequence` | `BioSequence{A}` | Typed biological sequence. |
| `identifier` | `String` | Stable id. |
| `name` | `String` | Short display name. |
| `description` | `String` | Longer description. |
| `annotations` | `Dict{Symbol,Any}` | Record-level annotations. |
| `letter_annotations` | `Dict{Symbol,Any}` | Per-position annotations such as quality scores. |

**Invariant:** String or vector values in `letter_annotations` must match `length(sequence)`.

---

## 4. Constructors

### `SeqRecord(sequence; kwargs...)`

```julia
SeqRecord(
    sequence::BioSequence{A};
    identifier="",
    description=identifier,
    metadata=Dict{Symbol,Any}(),
    prov_ctx=nothing
) where {A <: BioAlphabet}
```

**Behavior:**

- copies metadata into `Dict{Symbol,Any}`;
- ensures a provenance id;
- constructs `SeqRecord{A}`;
- records provenance when a context is active.

**Example:**

```julia
record = SeqRecord(DNASeq("ACGT"); identifier="seq1")
```

### `FastqRecord(sequence, quality; kwargs...)`

```julia
FastqRecord(
    sequence::BioSequence{A},
    quality::String;
    identifier="",
    description=identifier,
    metadata=Dict{Symbol,Any}(),
    prov_ctx=nothing
) where {A <: BioAlphabet}
```

**Behavior:**

- validates sequence and quality lengths;
- stores quality as raw string;
- ensures metadata provenance id;
- records provenance when active.

**Errors:** Throws `ArgumentError` if sequence and quality lengths differ.

**Example:**

```julia
fq = FastqRecord(DNASeq("ACGT"), "IIII"; identifier="read1")
```

### Backward-compatible FASTQ constructor

```julia
FastqRecord(identifier::String, description::String, sequence::BioSequence{A}, quality::String)
```

**Description:** Positional constructor used by older quality and FASTQ code paths.

### `SeqRecordLite(sequence::BioSequence; kwargs...)`

```julia
SeqRecordLite(
    sequence::BioSequence{A};
    identifier="",
    name=identifier,
    description=name,
    annotations=Dict{Symbol,Any}(),
    letter_annotations=Dict{Symbol,Any}(),
    prov_ctx=nothing
)
```

**Behavior:**

- validates letter annotation lengths;
- copies annotations;
- ensures provenance id in `annotations`;
- records provenance when active.

### `SeqRecordLite(sequence::AbstractString; kwargs...)`

```julia
SeqRecordLite(sequence::AbstractString; kwargs...)
```

**Description:** Infers the alphabet by checking DNA, then RNA, then amino-acid validity. Constructs a typed `BioSequence` and delegates to the typed constructor.

**Errors:** Throws `ArgumentError` when no alphabet can be inferred.

---

## 5. Utility Methods

### Length and conversion

```julia
length(record::Union{SeqRecord, FastqRecord, SeqRecordLite})
String(record::SeqRecordLite)
convert(String, record::SeqRecordLite)
```

**Description:** Record length delegates to sequence length. `SeqRecordLite` converts to its sequence string.

### Equality

```julia
==(left::SeqRecord{A}, right::SeqRecord{A})
==(left::FastqRecord{A}, right::FastqRecord{A})
==(left::SeqRecordLite{A}, right::SeqRecordLite{A})
```

**Description:** Records compare sequence and descriptive fields. Metadata comparison strips provenance keys so provenance stamping does not change biological equality.

### Display

```julia
show(record::SeqRecord)
show(record::FastqRecord)
show(record::SeqRecordLite)
```

**Example display:**

```text
SeqRecord{DNAAlphabet}(seq1, 4 bp, provenance=id=...)
FastqRecord{DNAAlphabet}(read1, 4 bp, provenance=id=...)
SeqRecordLite{AminoAcidAlphabet}(prot1, 120 bp, provenance=id=...)
```

---

## 6. Quality Scores

### `quality_scores`

```julia
quality_scores(record::FastqRecord; offset=33) -> Vector{Int}
```

**Description:** Converts the raw ASCII quality string to integer Phred scores.

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `offset` | `33` | ASCII offset for Phred encoding. |

**Example:**

```julia
fq = FastqRecord(DNASeq("ACGT"), "IIII")
quality_scores(fq)  # [40, 40, 40, 40]
```

---

## 7. Quick Reference

| API | Purpose |
|---|---|
| `SeqRecord{A}` | Annotated typed sequence record. |
| `FastqRecord{A}` | Typed sequence plus raw FASTQ quality string. |
| `SeqRecordLite{A}` | Lightweight annotated record with letter annotations. |
| `SeqRecord(...)` | Construct annotated sequence record. |
| `FastqRecord(...)` | Construct FASTQ record with length validation. |
| `SeqRecordLite(...)` | Construct lightweight record, optionally inferring alphabet from string. |
| `length(record)` | Sequence length. |
| `String(record::SeqRecordLite)` | Sequence string. |
| `quality_scores(record)` | Convert Phred+33 quality to integers. |

---

## Complete Usage Example

```julia
dna = DNASeq("ACGT")

rec = SeqRecord(dna; identifier="seq1", description="example DNA")
fq = FastqRecord(dna, "IIII"; identifier="read1")
lite = SeqRecordLite(
    dna;
    identifier="seq1",
    letter_annotations=Dict(:quality => "IIII"))

length(rec)
quality_scores(fq)
String(lite)
```
