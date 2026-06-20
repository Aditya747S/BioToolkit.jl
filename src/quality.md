# `quality.jl` - Sequencing Quality Utilities

## Overview

`quality.jl` provides FASTQ quality-score conversion, quality filtering, adapter trimming, low-quality trimming, and simple sequencing-record processing pipelines.

### Purpose

FASTQ workflows require consistent handling of Phred quality strings, quality thresholds, adapter removal, and minimum-length filtering. This file centralizes those operations for `FastqRecord` and `SeqRecordLite` while preserving provenance support.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Phred offset is explicit** | `offset=33` is the default, but callers can use other encodings. |
| **FASTQ and lightweight records both supported** | `FastqRecord` stores `quality`; `SeqRecordLite` stores quality in `letter_annotations`. |
| **Filtering returns booleans or records** | `quality_filter` answers pass/fail, while trimming/pipeline functions return new records or `nothing`. |
| **Adapter trimming supports both ends** | `_adapter_trim_span` can search from `:three_prime` or `:five_prime`. |
| **Provenance-aware helpers** | Public functions route through `active_provenance_context`. |

---

## 1. Phred Conversion

### `phred_scores`

```julia
phred_scores(quality::String; offset=33)
phred_scores(bytes::AbstractVector{UInt8}; offset=33)
phred_scores(record::FastqRecord; offset=33)
phred_scores(record::SeqRecordLite; quality_key=:quality, offset=33)
```

**Description:** Converts ASCII quality characters to integer Phred scores by subtracting `offset`.

**Example:**

```julia
phred_scores("IIII")  # [40, 40, 40, 40]
```

### `phred_string`

```julia
phred_string(scores::AbstractVector{<:Integer}; offset=33)
```

**Description:** Converts integer Phred scores back into an ASCII quality string.

---

## 2. Quality Summaries

### `mean_quality`

```julia
mean_quality(scores::AbstractVector{<:Integer}; offset=33)
mean_quality(quality::String; offset=33)
mean_quality(record::FastqRecord; offset=33)
mean_quality(record::SeqRecordLite; quality_key=:quality, offset=33)
```

**Description:** Computes mean Phred quality. Empty score vectors return `0.0`.

---

## 3. Filtering

### `quality_filter`

```julia
quality_filter(record; min_mean_quality=..., min_base_quality=..., max_low_quality_fraction=..., offset=33)
```

**Description:** Returns `true` when a record satisfies mean-quality and low-quality-base thresholds.

**Checks:**

- mean quality must be at least `min_mean_quality`;
- bases below `min_base_quality` are counted;
- low-quality fraction must be at most `max_low_quality_fraction`.

---

## 4. Adapter Trimming

### `_adapter_trim_span`

```julia
_adapter_trim_span(sequence, adapter; min_overlap=8, max_mismatches=1, from_end=:three_prime)
```

**Kind:** Internal helper

**Description:** Finds a trim span where an adapter overlaps either the three-prime or five-prime end. Mismatch tolerance is controlled by `max_mismatches`.

### `adapter_trim`

```julia
adapter_trim(record; adapter, min_overlap=8, max_mismatches=1, from_end=:three_prime)
```

**Description:** Trims adapter sequence from `FastqRecord` or `SeqRecordLite`, preserving quality/letter annotations consistently with the trimmed sequence.

---

## 5. Low-Quality Trimming

### `trim_low_quality`

```julia
trim_low_quality(record::FastqRecord; window=4, threshold=20, offset=33)
trim_low_quality(record::SeqRecordLite; quality_key=:quality, window=4, threshold=20, offset=33)
```

**Description:** Trims records using a sliding quality window. The helper `_trim_window` identifies retained sequence bounds.

---

## 6. Pipeline Helpers

### `process_sequencing_record`

```julia
process_sequencing_record(record; kwargs...) -> Union{record,Nothing}
```

**Description:** Applies adapter trimming, low-quality trimming, minimum-length filtering, and quality filtering to one record.

### `sequencing_pipeline`

```julia
sequencing_pipeline(records; kwargs...) -> Vector
```

**Description:** Processes a collection of sequencing records and returns records that pass all filters.

---

## Quick Reference

| API | Purpose |
|---|---|
| `phred_scores` | Quality string/record to integer scores. |
| `phred_string` | Integer scores to quality string. |
| `mean_quality` | Mean Phred quality. |
| `quality_filter` | Pass/fail quality thresholds. |
| `adapter_trim` | Trim adapter sequence. |
| `trim_low_quality` | Sliding-window quality trimming. |
| `process_sequencing_record` | Process one record. |
| `sequencing_pipeline` | Process a record collection. |

---

## Complete Usage Example

```julia
record = FastqRecord(DNASeq("ACGTACGT"), "IIII####"; identifier="read1")

scores = phred_scores(record)
mean_quality(record)

passes = quality_filter(record;
    min_mean_quality=20,
    min_base_quality=10,
    max_low_quality_fraction=0.5)

trimmed = trim_low_quality(record; window=4, threshold=20)
processed = process_sequencing_record(record; min_length=4)
```
