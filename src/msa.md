# msa.jl

## Purpose
This file defines the multiple sequence alignment data model used throughout BioToolkit. It provides the canonical container for aligned sequences and the metadata needed to reason about columns, conservation, and sequence-level annotations.

## Main structs
- AbstractMultipleSequenceAlignment: abstract parent for aligned sequence containers.
- MultipleSequenceAlignment: concrete alignment container with records, annotations, and column annotations.

## Public functions and methods
- Constructors for building a MultipleSequenceAlignment from sequence records.
- Record coercion helpers for converting different sequence-record types into alignment-friendly records.
- Column annotation validation and slicing helpers.
- Symbol-counting utilities for alignment columns.
- Conservation-related helpers that classify symbols across an aligned column.

## What the module does
The module stores aligned records in a way that preserves both the sequences and the metadata attached to the alignment. It lets downstream code inspect columns, count symbols, slice alignment regions, and keep annotations synchronized when the alignment is trimmed or reorganized.

## How the structs work together
The abstract alignment type gives BioToolkit a common API for aligned sequence sets. MultipleSequenceAlignment is the actual container users work with; it stores the aligned records and the associated annotations, including per-column metadata that is often important in comparative genomics or motif analysis.

## Typical usage
1. Create or load aligned sequence records.
2. Convert them into a MultipleSequenceAlignment.
3. Inspect or slice columns while preserving annotations.
4. Use symbol counting or conservation helpers to summarize columns.
5. Pass the alignment into downstream analysis or plotting routines.

## Important implementation details
- The module validates column annotations so the metadata stays aligned with the sequence length.
- Symbol counting helpers support downstream conservation summaries.
- The container is designed to work as a data model rather than a specific algorithm, so it can be reused across alignment analysis workflows.

## Threading notes
- `alignment_symbol_line()` now defaults to threaded column-wise execution when multiple threads are available.
- `consensus_sequence()` also defaults to threaded column-wise execution, with per-column counters allocated locally to keep the result deterministic.

## Why this file matters
MSA is a central data representation in comparative genomics and molecular evolution. This file gives BioToolkit a consistent place to store and manipulate aligned sequences before any downstream analysis or visualization.
