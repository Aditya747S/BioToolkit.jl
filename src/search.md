# search.jl

## Purpose
This file implements local sequence search and BLAST-style parsing. It gives BioToolkit a small k-mer index for internal searches, plus parsers and wrappers for standard BLAST outputs.

## Main types
- `KmerIndex` stores the packed k-mer lookup table, target names, and target sequences.
- `HighScoringPair` stores a local alignment hit in a BLAST-like format.
- `BlastXMLHSP`, `BlastXMLHit`, and `BlastXMLRecord` model BLAST XML output.
- `BlastTabularHSP`, `BlastTabularHit`, and `BlastTabularRecord` model BLAST tabular output.

## Public functions
- `build_index(targets; names, k)` builds a packed k-mer index.
- `local_search` is exported by the module for local matching workflows.
- `run_blast(program, query_seq, database; options, outfmt)` runs a local BLAST+ program and parses tabular hits.
- `qblast(program, database, query; options)` is a placeholder for a remote BLAST workflow.
- `parse_blast_xml`, `read_blast_xml`, `parse_blast_tabular`, and `read_blast_tabular` parse standard BLAST result formats.

## How it is used
The typical pattern is to build a k-mer index for a small target set, then use that index for quick local searches. When interoperability with existing tools matters, `run_blast` and the parsing helpers let the module ingest BLAST output rather than reimplement the whole search stack.

## Implementation notes
- k-mers are packed into `UInt64` values for efficient dictionary lookup.
- The BLAST wrapper writes temporary query/output files and shells out to the local BLAST+ binary.
- The XML and tabular parsers are intentionally lightweight and focused on the commonly used fields.

## Why it matters
Sequence search is one of the core daily tasks in bioinformatics. This file gives BioToolkit both a fast internal search primitive and a compatibility layer for BLAST-based workflows.
