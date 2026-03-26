# annotation.jl

## Purpose
This file defines the annotation model used by BioToolkit for gene features, sequence records, and feature slicing. It is the module that turns plain sequences into annotated biological records that can be queried, sliced, displayed, and converted between common feature formats.

## Main structs
- AbstractFeatureLocation: abstract parent for all feature-location objects.
- FeatureLocationLite: a simple interval with start, stop, strand, and partiality flags.
- CompoundFeatureLocation: a compound feature location made up of multiple parts such as joins or orders.
- SeqFeatureLite: a compact feature record with feature type, location, qualifiers, and identifier.
- AnnotatedSeqRecord: a mutable annotated sequence record with features, per-letter annotations, and metadata.
- SangerTrace: a chromatogram-style trace object with A, C, G, and T signal arrays plus qualities and annotations.

## Public functions and methods
- feature_spans(location or feature): returns the spans represented by a feature location.
- feature_bounds(location): returns the outer start and stop coordinates for a feature.
- feature_start and feature_stop: convenience accessors for the first and last coordinate.
- feature_strand: returns the strand direction encoded in a location.
- feature_identifier and feature_annotations: extract identifiers and qualifiers.
- feature_annotation(feature, key; default): read a single qualifier value.
- feature_contains: test whether a position falls inside a feature.
- feature_overlaps: test overlap between two features or locations.
- feature_extract: pull the subsequence for a location or feature.
- feature_slice: slice sequence records while preserving or transforming annotations.
- slice_feature_location and slice_feature: remap a feature into a sliced coordinate system.
- slice_annotated_record: slice a full annotated record.
- reverse_complement(record): reverse-complement an annotated record.
- annotate_gff_records and annotate_genbank_record(s): convert parsed format-specific records into SeqFeatureLite containers.
- annotate_variants: score and annotate variant records against gene features.
- parse_feature_location: parse GenBank/GFF-style location strings into typed locations.
- feature_table(record): convert annotated features into a tabular DataFrame.
- features_overlapping and select_features: query annotations by region or feature criteria.

## What the structs represent
FeatureLocationLite is the simplest case: a contiguous interval with orientation and optional partial ends. CompoundFeatureLocation handles complex feature expressions such as split exons or ordered intervals. SeqFeatureLite stores the feature type and qualifier dictionary in a form that is easy to inspect and manipulate. AnnotatedSeqRecord is the main high-level container, combining the raw sequence with annotations and features. SangerTrace is separate because it is a trace-format object rather than a gene-annotation object, but it lives here because it represents annotated sequence evidence.

## How the module is used
Typical use starts with parsing a record or converting a GenBank/GFF feature into SeqFeatureLite. From there, the user can:
1. Query feature coordinates with feature_start, feature_stop, or feature_bounds.
2. Slice records while keeping annotations synchronized.
3. Reverse-complement an annotated record when strand orientation changes.
4. Search features by overlap or containment.
5. Build a DataFrame with feature_table for reporting or downstream analysis.

## Internal design
The file contains a feature-location parser with caching, helper functions for splitting nested location expressions, and translation/variant-annotation helpers for consequence prediction. It also includes logic for handling GenBank and GFF-specific metadata so the package can normalize different annotation sources into the same data model.

## Why this file matters
This is the package's annotation bridge. It connects parsed file formats, biological feature semantics, and sequence slicing in a single place so downstream analysis code can work with one consistent record model.
