# record.jl

## Purpose
This file defines the canonical sequence-record abstractions used by BioToolkit. It normalizes simple FASTQ reads and lightweight annotated sequence records into a shared representation that other modules can read, write, and display consistently.

## Main structs
- FastqRecord: a plain FASTQ record with identifier, description, sequence, and quality string.
- SeqRecordLite: a mutable sequence record with identifier, name, description, annotations, and per-letter annotations.

## Public constructors and functions
- SeqRecordLite(sequence; identifier, name, description, annotations, letter_annotations): build a typed record from a raw sequence.
- SeqRecordLite(record::FastqRecord; annotations): convert a FASTQ record into a SeqRecordLite and preserve quality as a letter annotation.
- write_fastq(path, records; quality_key): export FASTQ records or SeqRecordLite objects to disk.
- Base.length(record::SeqRecordLite): return the sequence length.
- show methods for FastqRecord and SeqRecordLite: give readable REPL summaries.

## What the module does
The module provides a small but important layer of normalization. FASTQ data is converted into a lightweight object that can be annotated, while SeqRecordLite can store additional metadata without forcing a heavier sequence type. This makes downstream code simpler because it only has to understand a single record interface.

## How the structs work together
FastqRecord is the raw read format. SeqRecordLite is the more flexible in-memory record that can carry annotations and letter-level metadata such as qualities. When a FASTQ record is loaded into a SeqRecordLite, its quality string is preserved in letter_annotations so it can be written back out later.

## Typical usage
1. Read or construct FastqRecord objects when working with raw sequencing reads.
2. Convert them to SeqRecordLite if you need annotations or additional metadata.
3. Access length and display methods for quick inspection.
4. Call write_fastq when you want to export reads back to FASTQ.

## Important implementation details
- _fastq_components extracts the tuple of FASTQ fields from either record type.
- write_fastq enforces that sequence and quality lengths match.
- SeqRecordLite stores annotations and letter_annotations as dictionaries, which makes it flexible enough for many bioinformatics workflows.

## Why this file matters
This module is the package's basic sequence-record abstraction layer. It gives the rest of BioToolkit a stable object to use for raw reads and lightweight annotated sequences.
