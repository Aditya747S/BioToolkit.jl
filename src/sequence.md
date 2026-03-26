# sequence.jl

## Purpose
This file provides the low-level DNA and nucleotide sequence primitives used throughout BioToolkit. It is built around fast lookup tables for sequence encoding, complements, translation, and other operations that need to be fast enough for large-scale sequence analysis.

## Main structs
- FastaIndexRecord: a compact metadata record for indexed FASTA files.

## Public functions and helpers
- sequence parsing and encoding helpers for DNA and IUPAC symbols.
- complement and reverse-complement logic for nucleotide strings.
- codon translation helpers based on the embedded codon table.
- molecular weight and residue lookup helpers for nucleic acid calculations.
- FASTA indexing support through FastaIndexRecord.

## What the module does
The module provides the basic operations needed to work with DNA strings efficiently. It stores the alphabet encoding in lookup tables so conversion from characters to codes is fast, and it also provides codon and complement tables for downstream annotation and alignment modules.

## How the struct is used
FastaIndexRecord is used when a FASTA file needs compact random-access metadata. The rest of the module is table-driven rather than object-heavy, which makes it suitable for tight loops and repeated transformations.

## Typical usage
1. Encode or normalize a DNA string into the package's internal nucleotide representation.
2. Complement or reverse-complement sequences when working with strands.
3. Translate codons into amino-acid sequences for protein-oriented workflows.
4. Use the molecular-weight and base lookup tables when you need fast physical calculations.

## Important implementation details
- _DNA_CODE stores a character-to-code mapping for DNA and IUPAC symbols.
- _DNA_COMPLEMENT and _CODON_TABLE provide precomputed biology-specific lookups.
- The module is deliberately optimized around constant-time table access instead of repeated parsing or string conversion.
- CUDA-backed sequence helpers are now opt-in through a `use_cuda=false` keyword, so the CPU path remains the default.
- When `use_cuda=true`, the public wrappers promote CPU inputs to CUDA arrays before calling the lazy GPU kernels, so callers do not need to pre-wrap inputs themselves.

## Why this file matters
This module is one of the most performance-sensitive foundations in BioToolkit. It keeps the most common sequence operations fast and predictable so higher-level code can build on it without paying repeated conversion costs.
