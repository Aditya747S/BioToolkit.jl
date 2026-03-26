# bam.jl

## Purpose
This file defines the BAM data model and most of the low-level binary encoding and decoding logic for alignment files. It is responsible for representing BAM headers, records, indices, and the conversion between packed BAM data and Julia objects.

## Main structs
- BamReference: a reference sequence name and length entry from the BAM header.
- BamHeader: the textual header plus its list of references.
- BamCigarOp: a single CIGAR operation and its length.
- BamRecord: one alignment record containing read name, flags, reference positions, CIGAR operations, mate information, sequence, qualities, and tags.
- BamFile: a simple in-memory container holding a header and a list of records.
- BamChunk: a compressed file chunk used for random-access indexing.
- BamIndex: the binning and linear index tables used for region queries.

## Public constructors and Base methods
- BamReference(name, length): normalize a reference entry.
- BamHeader(references; text=""): create or infer a header string from references.
- BamHeader(text, references): construct a header from explicit text and references.
- BamCigarOp(length, op): create a CIGAR operation.
- BamRecord(...): construct a fully typed alignment record.
- BamFile(records; header=nothing): wrap a list of records in a file container and infer the header if needed.
- == methods for BamReference, BamCigarOp, BamRecord, BamHeader, and BamFile.
- show methods for readable summaries of the BAM-related structs.

## Important internal helpers
- _bam_header_text: generate a minimal SAM-style header text block.
- _decode_bam_sequence and _encode_bam_sequence: translate between BAM nucleotide encoding and plain strings.
- _decode_bam_quality and _encode_bam_quality: convert quality data.
- _decode_bam_cigar and _encode_bam_cigar: convert CIGAR binary codes to and from BamCigarOp.
- _parse_bam_tags and _write_bam_tags: handle auxiliary tag fields.
- _bam_reference_name_map, _reference_index, and _reference_name: resolve reference IDs and names.
- _bam_reference_span and _bam_query_span: compute read and alignment spans.
- _bam_overlaps: test overlap between a record and a genomic interval.
- _bam_record_from_stream: decode a BAM record from a BGZF stream.
- _bam_reg2bins and _bam_region_chunks: compute index bins and chunks for regional queries.
- _infer_bam_header and _bam_infer_reference_lengths: infer header content when it is missing.
- _write_bam_record and _write_bam_index: serialize BAM data and its index.
- _read_bam_index and _read_bam_header: load BAM metadata from disk.
- _bam_scan_region and _bam_region_from_index: read alignments restricted to a region.

## Public I/O functions
- read_bam(path): read a BAM file into a BamFile object.
- read_bam(path, region): read only alignments overlapping a genomic interval.
- write_bam(path, records; header, write_index): write records to disk and optionally generate an index.

## Typical usage
1. Read a BAM file with read_bam(path).
2. Inspect the returned BamFile, its header, and its records.
3. Query a region with read_bam(path, region) when you only need a genomic slice.
4. Construct BamRecord objects manually when exporting data or creating synthetic alignments.
5. Use write_bam to save edited or generated records back to disk.

## Why this file matters
This is the foundational BAM layer for the package. It keeps the binary format details in one place so the rest of BioToolkit can work with high-level alignment records instead of byte-level BAM structures.
