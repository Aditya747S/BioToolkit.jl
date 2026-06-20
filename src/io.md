# `io.jl` - Biological Format I/O

## Overview

`io.jl` implements BioToolkit's core biological file parsers, writers, and ingestion helpers. It covers sequence formats, variant and interval formats, GenBank/EMBL-like records, Arrow ingestion, and selected instrument formats.

### Purpose

Bioinformatics code spends a large amount of time moving data between files and typed in-memory objects. `io.jl` provides high-performance, provenance-aware readers and writers for common formats:

- FASTA and FASTQ sequence records;
- VCF documents and variant records;
- BED and GFF interval records;
- GenBank and EMBL sequence annotations;
- Arrow table ingestion for larger tabular genomics workflows;
- ABIF trace reads;
- SwissProt placeholder/compatibility entry points.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Typed alphabet dispatch** | FASTA/FASTQ readers produce `SeqRecord{A}` or `FastqRecord{A}` using `BioSequence{A}`. |
| **Byte-level parsing where useful** | Sequence and trace formats avoid unnecessary high-level parsing overhead. |
| **Path and IO overloads** | Most readers support file paths and existing `IO` streams. |
| **Provenance hashes for path reads/writes** | Path-based reads hash input bytes; writes hash output files after writing. |
| **Optional gzip VCF support** | `.vcf.gz` support loads `CodecZlib` only when needed. |
| **Chunked Arrow ingestion** | Large VCF/BED/GFF/GenBank inputs can be streamed into Arrow chunks. |

---

## Table of Contents

1. [FASTA I/O](#1-fasta-io)
2. [FASTQ I/O](#2-fastq-io)
3. [VCF I/O](#3-vcf-io)
4. [BED and GFF I/O](#4-bed-and-gff-io)
5. [GenBank I/O](#5-genbank-io)
6. [Arrow Ingestion](#6-arrow-ingestion)
7. [EMBL, SwissProt, and ABIF](#7-embl-swissprot-and-abif)
8. [Arrow Table Helpers](#8-arrow-table-helpers)
9. [Quick Reference](#9-quick-reference)

---

## 1. FASTA I/O

### `read_fasta`

```julia
read_fasta(path::String; alphabet=DNAAlphabet, prov_ctx=nothing)
read_fasta(io::IO; alphabet=DNAAlphabet, prov_ctx=nothing)
```

**Description:** Reads FASTA records into `Vector{SeqRecord{A}}`, where `A` is the requested alphabet type.

**Behavior:**

- path reads load raw bytes and compute SHA-256 provenance hash;
- `>` starts a new record;
- blank lines are ignored;
- sequence lines are uppercased;
- sequence data before the first header throws `ArgumentError`;
- records are built with `BioSequence{A}(bytes; validate=false)`;
- provenance nodes are registered for the file and each record when active.

**Example:**

```julia
records = read_fasta("genome.fa"; alphabet=DNAAlphabet)
records[1].identifier
records[1].sequence
```

### `write_fasta`

```julia
write_fasta(path::String, records; prov_ctx=nothing) -> path
```

**Description:** Writes FASTA records with `>` headers and 60-character sequence lines.

**Expected record shape:** Each record must expose `identifier` and `sequence.data`.

---

## 2. FASTQ I/O

### `read_fastq`

```julia
read_fastq(path::String; alphabet=DNAAlphabet, prov_ctx=nothing)
read_fastq(io::IO; alphabet=DNAAlphabet, prov_ctx=nothing)
```

**Description:** Reads FASTQ records into `Vector{FastqRecord{A}}`.

**Validation:**

- first line of each record must start with `@`;
- third line must start with `+`;
- sequence construction validates against `alphabet`;
- `FastqRecord` enforces sequence/quality length equality.

**Identifier behavior:** `_fastq_identifier` takes the first whitespace-delimited token from the header as `identifier`, while the full header becomes `description`.

### `write_fastq`

```julia
write_fastq(path::String, records; prov_ctx=nothing) -> path
```

**Description:** Writes FASTQ records using `_fastq_components(record)` to support `FastqRecord` and compatible `SeqRecordLite` values.

---

## 3. VCF I/O

### `parse_vcf_record`

```julia
parse_vcf_record(line::AbstractString; prov_ctx=nothing)
```

**Description:** Parses one non-header VCF line into a `VariantTextRecord`. Returns `nothing` for records with too few fields when appropriate; malformed records throw `ArgumentError`.

### `read_vcf_document`

```julia
read_vcf_document(input::String; prov_ctx=nothing)
read_vcf_document(io::IO; prov_ctx=nothing)
```

**Description:** Reads a full VCF document into `VcfDocument`, preserving header meta-lines, column names, sample names, and variant records.

**Path behavior:** `.vcf.gz` inputs are supported through optional `CodecZlib`.

### `read_vcf`

```julia
read_vcf(input::String)
read_vcf(io::IO)
```

**Description:** Convenience reader returning VCF records/document data through the package's VCF parsing path.

### `write_vcf_document`

```julia
write_vcf_document(output::String, doc::VcfDocument)
write_vcf_document(io::IO, doc::VcfDocument)
```

**Description:** Writes a full VCF document with header and records.

### `write_vcf`

```julia
write_vcf(output::String, records; header=nothing)
write_vcf(io::IO, records; header=nothing)
write_vcf(output::String, doc::VcfDocument)
write_vcf(io::IO, doc::VcfDocument)
```

**Description:** Writes variant records or a `VcfDocument` to VCF text.

**Important behavior:** Sample count mismatches in records raise `DimensionMismatch`.

---

## 4. BED and GFF I/O

### `BedRecord`

```julia
struct BedRecord
    chrom
    start
    stop
    ...
end
```

**Description:** BED interval record. BED uses 0-based starts and half-open stops at the file-format level; downstream conversion to `GenomicInterval` adjusts start by +1.

### `GffRecord`

```julia
struct GffRecord
    chrom
    source
    feature
    start
    stop
    score
    strand
    phase
    attributes
    attribute_map
end
```

**Description:** GFF3-like feature record with parsed attributes.

### BED API

```julia
parse_bed_record(line)
read_bed(input)
write_bed(output, records)
```

**Behavior:**

- comment lines beginning with `#` are skipped;
- malformed data include line-number context in reader errors;
- writers emit at least chromosome, start, and stop.

### GFF API

```julia
parse_gff_record(line)
read_gff(input)
write_gff(output, records)
```

**Behavior:**

- parses attributes into a `Dict{String,Vector{String}}`;
- preserves explicit attribute text when present;
- renders attribute maps when raw attributes are empty.

---

## 5. GenBank I/O

### `GenBankFeature`

```julia
struct GenBankFeature
    key::String
    location::String
    qualifiers::Dict{String,Vector{String}}
    parsed_location
end
```

**Description:** One GenBank feature table entry with raw and parsed location information.

### `GenBankRecord`

```julia
struct GenBankRecord
    locus
    locus_line
    definition
    accession
    version
    keywords
    source
    organism
    features
    sequence
    comment
    metadata
end
```

**Description:** Parsed GenBank sequence record with annotation and sequence data.

### `parse_genbank_record`

```julia
parse_genbank_record(lines::AbstractVector{<:String})
```

**Description:** Parses one GenBank record from its text lines. Handles locus, definition, accession, version, keywords, source, organism, comment, features, qualifiers, and origin sequence.

### `read_genbank`

```julia
read_genbank(input_path::String; prov_ctx=nothing)
read_genbank(io::IO; prov_ctx=nothing)
```

**Description:** Reads GenBank records separated by `//`.

### `write_genbank`

```julia
write_genbank(output_path::String, records; prov_ctx=nothing)
write_genbank(io::IO, records; prov_ctx=nothing)
```

**Description:** Writes GenBank records, wrapping text sections and rendering features and sequence origin blocks.

---

## 6. Arrow Ingestion

### `ingest_genbank`

```julia
ingest_genbank(input_path, output_path; chunk_size=100)
```

**Description:** Streams GenBank records into an Arrow file in chunks.

### `ingest_vcf`

```julia
ingest_vcf(input_path, output_path; chunk_size=10_000)
```

**Description:** Streams VCF rows into Arrow columns, including sample-count and sample text columns.

### `ingest_bed`

```julia
ingest_bed(input_path, output_path; chunk_size=10_000)
```

**Description:** Streams BED rows into Arrow columns.

### `ingest_gff`

```julia
ingest_gff(input_path, output_path; chunk_size=10_000)
```

**Description:** Streams GFF rows into Arrow columns.

---

## 7. EMBL, SwissProt, and ABIF

### `read_embl`

```julia
read_embl(filepath::String)
read_embl(io::IO)
```

**Description:** Reads EMBL records and converts them into lightweight sequence records with feature annotations.

### `read_swissprot`

```julia
read_swissprot(filepath::String)
```

**Description:** SwissProt reader entry point. See source for current implementation details and limitations.

### `read_abif`

```julia
read_abif(filepath::String)
read_abif(io::IO)
```

**Description:** Reads ABIF trace files, validates the `ABIF` magic, parses directory entries, and extracts trace channels where present.

---

## 8. Arrow Table Helpers

### `load_arrow_table`

```julia
load_arrow_table(path::String)
```

**Description:** Loads an Arrow table from disk.

### `write_arrow_table`

```julia
write_arrow_table(output_path::String, table)
```

**Description:** Writes a table to Arrow format.

---

## 9. Quick Reference

| API | Purpose |
|---|---|
| `read_fasta` / `write_fasta` | FASTA records. |
| `read_fastq` / `write_fastq` | FASTQ records. |
| `parse_vcf_record` | Parse one VCF data line. |
| `read_vcf_document` / `write_vcf_document` | Full VCF documents. |
| `read_vcf` / `write_vcf` | VCF convenience API. |
| `parse_bed_record`, `read_bed`, `write_bed` | BED records. |
| `parse_gff_record`, `read_gff`, `write_gff` | GFF records. |
| `read_genbank`, `write_genbank` | GenBank records. |
| `ingest_vcf`, `ingest_bed`, `ingest_gff`, `ingest_genbank` | Chunked Arrow ingestion. |
| `read_embl` | EMBL records. |
| `read_abif` | ABIF trace files. |
| `load_arrow_table`, `write_arrow_table` | Arrow table helpers. |

---

## Complete Usage Example

```julia
records = read_fasta("input.fa"; alphabet=DNAAlphabet)
write_fasta("copy.fa", records)

fastq = read_fastq("reads.fastq"; alphabet=DNAAlphabet)
scores = quality_scores(first(fastq))

vcf_doc = read_vcf_document("variants.vcf")
write_vcf_document("variants.copy.vcf", vcf_doc)

bed = read_bed("regions.bed")
gff = read_gff("annotation.gff3")

gb = read_genbank("records.gb")
write_genbank("records.copy.gb", gb)
```
