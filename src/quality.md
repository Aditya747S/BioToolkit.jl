# quality.jl

## Purpose
This file handles FASTQ quality-score conversion, read trimming, and basic quality control. It is the small preprocessing layer that sits between raw sequencing input and downstream analysis.

## Main functions
- `phred_score(byte; offset=33)` converts a single encoded quality byte to a Phred score.
- `phred_scores(quality; offset=33)` converts a string or byte vector into numeric scores.
- `phred_scores(record; offset=33)` supports both `FastqRecord` and `SeqRecordLite` inputs.
- `phred_string(scores; offset=33)` converts numeric scores back into an ASCII quality string.
- `mean_quality(scores; offset=33)` computes the average Phred score.
- `quality_filter(record; min_mean_quality, min_base_quality, max_low_quality_fraction, quality_key, offset)` applies a read-level filter.
- `trim_low_quality(record; window, threshold, quality_key, offset)` trims low-quality ends from a record.
- `adapter_trim(record; adapter, min_overlap, max_mismatches, from_end, quality_key)` removes adapter contamination from either end.
- `process_sequencing_record(record; ...)` combines adapter trimming, end trimming, length filtering, and quality filtering.
- `sequencing_pipeline(records; ...)` applies the same processing logic to a vector of records.

## Internal helpers
- `_quality_scores` normalizes the record-specific quality extraction path.
- `_trim_window` identifies the longest retained high-quality span.
- `_adapter_trim_span` finds the adapter overlap to trim from either the 3' or 5' end.

## How it is used
The normal usage pattern is to take a raw FASTQ-like record, convert or inspect the quality encoding with `phred_scores` or `mean_quality`, then remove low-quality sequence with `trim_low_quality` or `process_sequencing_record`. When many reads need the same treatment, `sequencing_pipeline` applies the same rules to a batch.

The functions are written to work with both `FastqRecord` and `SeqRecordLite`, so the module can sit on top of either the FASTQ-specific path or the more general record abstraction.

## Implementation notes
- Quality values are assumed to use a configurable ASCII offset, with 33 as the default.
- Adapter trimming searches for a best overlap within the allowed mismatch budget and only trims when a valid match is found.
- `quality_filter` rejects empty reads and can enforce mean-quality, base-quality, and low-quality-fraction constraints together.

## Why it matters
Quality handling is one of the first gates in a sequencing workflow. This file centralizes the common cleanup logic so later alignment, variant, or assembly steps do not need to duplicate it.
