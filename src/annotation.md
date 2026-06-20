# `annotation.jl` - Sequence Annotation Utilities

## Overview

`annotation.jl` provides feature-location types, annotated sequence records, feature slicing, GenBank/GFF conversion, feature tables, and simple variant consequence annotation.

### Purpose

This document is a source-matched reference for the public API in `annotation.jl`. It groups exported types and functions by workflow stage and explains the biological or statistical role of each entry.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Workflow-oriented grouping** | Functions are organized by how users run the analysis. |
| **Typed domain state** | Reusable records and result objects are represented explicitly. |
| **Scalable summaries** | APIs return compact matrices/tables for large biological datasets. |
| **Interoperable outputs** | Results can feed plotting, enrichment, reporting, and downstream modeling modules. |
| **Explicit thresholds and controls** | Filtering, QC, and model functions expose important parameters as keywords. |

---

## 1. Feature Model

Annotation objects describe single intervals, compound joins/orders, feature qualifiers, and sequence-level metadata.

| API | Description |
|---|---|
| `AbstractFeatureLocation` | Common supertype for simple and compound feature locations. |
| `FeatureLocationLite` | Simple interval with start, stop, strand, and partial-end flags. |
| `CompoundFeatureLocation` | Joined/ordered collection of feature-location parts. |
| `SeqFeatureLite` | Feature type, location, qualifiers, and feature id. |
| `AnnotatedSeqRecord` | Sequence plus id/name/description, features, annotations, and letter annotations. |
| `SangerTrace` | Container for Sanger trace/channel data and base calls. |

## 2. Location Queries

These helpers compute bounds, spans, strand, containment, and overlap relationships.

| API | Description |
|---|---|
| `parse_feature_location` | Parses GenBank-style location strings, including complements and joins. |
| `feature_spans` | Returns one or more concrete spans for a location or feature. |
| `feature_start` | Returns the minimum start coordinate. |
| `feature_stop` | Returns the maximum stop coordinate. |
| `feature_bounds` | Returns start/stop bounds for a location. |
| `feature_length` | Computes total feature length across parts. |
| `feature_strand` | Returns feature strand, respecting compound locations. |
| `feature_contains` | Tests whether a position falls inside a feature/location. |
| `feature_overlaps` | Tests whether two features/locations overlap. |

## 3. Slicing and Extraction

Slicing functions preserve or transform feature coordinates and annotations.

| API | Description |
|---|---|
| `feature_sequence` | Extracts sequence covered by a feature/location. |
| `feature_extract` | Alias-style extraction for records or sequences. |
| `feature_slice` | Slices `SeqRecordLite`, `AnnotatedSeqRecord`, or `GenBankRecord` over a coordinate range. |
| `slice_feature_location` | Transforms a location into slice-relative coordinates. |
| `slice_feature` | Transforms a feature into slice-relative coordinates. |
| `slice_annotated_record` | Slices sequence, features, annotations, and letter annotations together. |
| `reverse_complement` | Reverse-complements annotated records while transforming feature coordinates. |

## 4. Selection, Tables, and Format Conversion

Converters normalize external formats into lightweight BioToolkit annotations.

| API | Description |
|---|---|
| `feature_summary` | Returns a compact named-tuple summary for one feature. |
| `select_features` | Filters features by type, id, qualifier, strand, or region. |
| `features_at` | Returns features overlapping one position. |
| `features_overlapping` | Returns features overlapping a location/range. |
| `feature_table` | Builds a DataFrame-style feature table. |
| `feature_annotations` | Returns feature qualifier dictionary. |
| `feature_annotation` | Retrieves one qualifier with a default. |
| `annotate_gff_records` | Converts GFF records to `SeqFeatureLite`s. |
| `annotate_genbank_record` | Converts a GenBank record to an annotated record. |
| `annotate_genbank_records` | Batch GenBank annotation conversion. |
| `annotate_variants` | Annotates variants against gene features and optional reference sequences. |

---

## Complete Usage Example

```julia
using BioToolkit

loc = parse_feature_location("complement(join(10..20,30..40))")
feature = SeqFeatureLite("CDS", loc; qualifiers=Dict("gene"=>["abc"]))
record = AnnotatedSeqRecord(DNASeq("ACGT"^20); features=[feature])
cds = feature_sequence(record, feature)
table = feature_table(record)
```

