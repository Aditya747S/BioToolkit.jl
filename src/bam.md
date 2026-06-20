# `bam.jl` - Pure-Julia BAM I/O

## Overview

`bam.jl` provides BAM headers, records, CIGAR operations, BGZF-backed streaming readers, materialized BAM files, BAM writing, and BAI-style index support.

### Purpose

This module lets BioToolkit read, create, scan, and write coordinate-aligned sequencing records without relying on htslib for the common BAM path.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Domain-specific containers** | Results use explicit structs where the workflow has stable fields or needs provenance metadata. |
| **Workflow APIs** | Functions expose complete analysis steps, not only internal numeric kernels. |
| **Table-compatible outputs** | Outputs are designed to work with Julia arrays, dictionaries, named tuples, and DataFrames. |
| **Pure-Julia core** | External tools are optional wrappers; core summaries remain usable inside Julia. |
| **Provenance hooks** | Many public functions accept provenance context keywords or return result containers with provenance records. |

---

## 1. Core Types

These types model BAM metadata, alignments, files, and index chunks.

| API | Description |
|---|---|
| `BamReference` | Reference sequence name and length from the BAM header. |
| `BamHeader` | Header text plus ordered `BamReference` records. |
| `BamCigarOp` | Single CIGAR operation with operation character and length. |
| `BamRecord` | One alignment record with query name, flag, reference, position, CIGAR, sequence, quality, mate fields, and tags. |
| `BamFile` | Materialized BAM container holding a header and records. |
| `BamChunk` | Virtual-offset span used by BAI bins. |
| `BamIndex` | In-memory BAM index with bins and linear offsets. |

## 2. Reading

Readers support full-file iteration and region-restricted access.

| API | Description |
|---|---|
| `BamReader` | Streaming reader over all records in a BAM file. |
| `BamRegionScanReader` | Region reader that scans records and filters overlaps when no index path is used. |
| `BamIndexedRegionReader` | Region reader that seeks through BAI chunks for indexed queries. |
| `BamEmptyReader` | Empty iterator returned for regions with no matching chunks. |
| `read_bam` | Reads a BAM file, optionally materializing records or restricting to a `GenomicInterval`. |

## 3. Writing and Indexing

Writers encode records, CIGARs, tags, qualities, and optional BAI indexes.

| API | Description |
|---|---|
| `write_bam` | Writes records to a BGZF-compressed BAM file with an inferred or supplied header. |
| `write_bam_index` | Writes a BAI-compatible index for a coordinate-sorted BAM. |
| `read_bam_index` | Loads BAI bins/chunks and linear indexes from disk. |

---

## Quick Reference

| Area | Main APIs |
|---|---|
| Core Types | `BamReference`, `BamHeader`, `BamCigarOp`, `BamRecord`, `BamFile`, `BamChunk`, ... |
| Reading | `BamReader`, `BamRegionScanReader`, `BamIndexedRegionReader`, `BamEmptyReader`, `read_bam` |
| Writing and Indexing | `write_bam`, `write_bam_index`, `read_bam_index` |

---

## Complete Usage Example

```julia
using BioToolkit

ref = BamReference("chr1", 1_000_000)
header = BamHeader([ref])
record = BamRecord("read1", DNASeq("ACGT"); refname="chr1", pos=0, cigar=[BamCigarOp(4, M)])
write_bam("example.bam", [record]; header=header)
bam = read_bam("example.bam")
```

