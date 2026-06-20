# `sequence.jl` - Nucleotide Sequence Operations

## Overview

`sequence.jl` provides BioToolkit's core DNA/RNA sequence operations: validation, FASTQ streaming, nucleotide counting, transcription, reverse complementation, GC analysis, melting temperature, molecular weight, codon usage, translation, k-mer analysis, FASTA indexing, primer/PCR utilities, skew analysis, dot matrices, and sequence complexity metrics.

### Purpose

This module is the main algorithm layer above `biotypes.jl`. It uses typed `BioSequence` inputs for alphabet safety and byte-indexed lookup tables for fast per-base operations.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Typed sequence dispatch** | DNA, RNA, and amino-acid operations can be separated by method signatures. |
| **256-element byte lookup tables** | Complement, codon, mass, and validation lookups are O(1). |
| **String overloads for compatibility** | Many functions still accept strings and convert to typed sequences internally. |
| **In-place APIs where useful** | Translation and reverse complement have buffer-writing variants. |
| **Optional CUDA paths** | Hamming, k-mer, and dot-matrix operations can use lazy CUDA helpers. |
| **FASTA `.fai` compatibility** | Index records support random access into FASTA files. |

---

## 1. FASTA Indexing and FASTQ Streaming

### `FastaIndexRecord`

```julia
struct FastaIndexRecord
    name::String
    sequence_length::Int
    byte_offset::Int64
    line_bases::Int
    line_bytes::Int
end
```

**Description:** Samtools-compatible `.fai` metadata for random access into FASTA files.

### `fasta_index`

```julia
fasta_index(path::String)
```

**Description:** Builds or reads index-style metadata for FASTA random access.

### `fetch_fasta_sequence`

```julia
fetch_fasta_sequence(path_or_bytes, index, start_pos, stop_pos)
```

**Description:** Fetches a subsequence from indexed FASTA data using byte offsets and line geometry.

### `FastqRecordStream` and `each_fastq_record`

```julia
each_fastq_record(path; sequence_type=DNASeq)
each_fastq_record(io; sequence_type=DNASeq)
```

**Description:** Streaming FASTQ iterator that yields records without loading the whole file at once.

---

## 2. Basic Nucleotide Operations

### `validate_dna`

```julia
validate_dna(sequence) -> Bool
```

**Description:** Checks whether a typed sequence or string is DNA-compatible.

### `count_nucleotides`

```julia
count_nucleotides(sequence) -> NamedTuple
```

**Description:** Counts A, C, G, T, and N-like bases.

### `transcribe_dna`

```julia
transcribe_dna(sequence)
```

**Description:** Converts DNA thymine to RNA uracil and returns an RNA sequence for typed DNA input.

### `reverse_complement`

```julia
reverse_complement(sequence)
reverse_complement(bytes)
reverse_complement!(buffer, bytes)
```

**Description:** Computes reverse complement for DNA/RNA/byte inputs. The in-place form writes into `buffer`.

---

## 3. Composition and Physical Properties

### `gc_content`

```julia
gc_content(sequence) -> Float64
```

**Description:** Fraction of G/C bases. Empty sequences return `0.0`. Unsupported symbols can throw `ArgumentError`.

### `melting_temp`

```julia
melting_temp(sequence; rna=false, method=:wallace)
```

**Description:** Estimates melting temperature using supported formulas such as the Wallace rule.

### `dna_molecular_weight`

```julia
dna_molecular_weight(sequence; stranded=:single)
```

**Description:** Computes DNA molecular weight using residue mass tables.

---

## 4. Codons and Translation

### `codon_usage`

```julia
codon_usage(sequence; frame=1, include_stop=true)
codon_usage(sequences; frame=1, include_stop=true)
```

**Description:** Counts codons in one DNA sequence or a vector of DNA sequences.

### `codon_usage_table`

```julia
codon_usage_table(sequence; frame=1, include_stop=false)
```

**Description:** Returns codon usage in table-like form.

### `relative_codon_adaptiveness`

```julia
relative_codon_adaptiveness(reference; frame=1, pseudocount=0.5)
```

**Description:** Builds reference codon weights for codon adaptation index calculations.

### `codon_adaptation_index` / `cai`

```julia
codon_adaptation_index(sequence; reference, frame=1, pseudocount=0.5)
cai(sequence; reference, kwargs...)
```

**Description:** Computes Sharp-Li-style codon adaptation index.

### `translate_dna`

```julia
translate_dna(sequence; stop_at_stop=false)
translate_dna(bytes; stop_at_stop=false)
translate_dna!(buffer, bytes; stop_at_stop=false)
```

**Description:** Translates DNA codons to amino-acid sequence. In-place forms require a sufficiently large buffer.

---

## 5. Distance, K-mers, and Dot Matrices

### `hamming_distance`

```julia
hamming_distance(left, right; use_cuda=false)
```

**Description:** Counts mismatches between equal-length sequences.

**Errors:** Throws `ArgumentError` if lengths differ.

### `kmer_frequency`

```julia
kmer_frequency(sequence, k; use_cuda=false)
```

**Description:** Counts k-mers in typed sequences, strings, or byte vectors. DNA has optimized paths for unambiguous and ambiguous input.

### `dotmatrix`

```julia
dotmatrix(seq1, seq2; use_cuda=false)
```

**Description:** Builds a match matrix between two sequences.

---

## 6. ORFs, PCR, and Guide Utilities

### `find_orfs`

```julia
find_orfs(sequence; min_aa=0)
```

**Description:** Finds translated open reading frames in DNA sequence.

### `protein_search`

```julia
protein_search(sequence, motif)
```

**Description:** Translates DNA ORFs and searches for amino-acid motif hits.

### `check_primer_dimers`

```julia
check_primer_dimers(fwd_seq, rev_seq; min_3prime_match=4)
```

**Description:** Detects potential 3-prime primer dimer matches.

### `pcr_in_silico`

```julia
pcr_in_silico(primer_fwd, primer_rev, genome_index; genome_path=nothing, max_product_length=50000)
```

**Description:** Searches genome sequences or indexed FASTA records for amplicons.

### `score_grna_efficiency`

```julia
score_grna_efficiency(guide; pam=DNASeq("NGG"))
```

**Description:** Computes a simple guide RNA efficiency score.

---

## 7. Skew and Complexity

### `gc_skew` and `minimum_skew`

```julia
gc_skew(sequence)
minimum_skew(sequence)
```

**Description:** Computes cumulative GC skew and positions of minimum skew.

### Frequency helpers

```julia
dinucleotide_frequency(sequence)
trinucleotide_frequency(sequence)
letter_frequency(sequence)
oligonucleotide_frequency(sequence; k=2, normalize=false)
```

**Description:** Counts sequence words or letters.

### Complexity helpers

```julia
sequence_complexity(sequence; k=3)
sequence_entropy(sequence; k=1)
sliding_window_gc_content(sequence; window=100, step=1)
```

**Description:** Computes linguistic complexity, entropy, and sliding-window GC content.

---

## Quick Reference

| Area | APIs |
|---|---|
| FASTA/FASTQ | `FastaIndexRecord`, `fasta_index`, `fetch_fasta_sequence`, `each_fastq_record` |
| Basic DNA/RNA | `validate_dna`, `count_nucleotides`, `transcribe_dna`, `reverse_complement` |
| Composition | `gc_content`, `melting_temp`, `dna_molecular_weight` |
| Codons | `codon_usage`, `codon_usage_table`, `relative_codon_adaptiveness`, `codon_adaptation_index`, `cai` |
| Translation | `translate_dna`, `translate_dna!` |
| Distance/k-mers | `hamming_distance`, `kmer_frequency`, `dotmatrix` |
| ORF/PCR | `find_orfs`, `protein_search`, `check_primer_dimers`, `pcr_in_silico` |
| CRISPR utility | `score_grna_efficiency` |
| Skew/complexity | `gc_skew`, `minimum_skew`, `sequence_complexity`, `sequence_entropy` |

---

## Complete Usage Example

```julia
dna = DNASeq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

count_nucleotides(dna)
gc_content(dna)
reverse_complement(dna)
translate_dna(dna)

codons = codon_usage(dna)
kmers = kmer_frequency(dna, 3)

orfs = find_orfs(dna; min_aa=10)
skew = gc_skew(dna)
```
