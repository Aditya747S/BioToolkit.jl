# `align.jl` - Pairwise and Codon-Aware Sequence Alignment

## Overview

`align.jl` implements global, local, affine-gap, matrix-scored, and codon-aware pairwise sequence alignment.

### Purpose

BioToolkit needs a native alignment layer that works directly with `BioSequence` values and byte vectors. This module provides Needleman-Wunsch global alignment, Smith-Waterman local alignment, Gotoh affine-gap alignment, named substitution matrices, custom matrices, and codon-space alignment.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Typed result object** | `PairwiseAlignmentResult{A}` stores aligned `BioSequence{A}` values plus score/identity. |
| **Linear and matrix scoring** | Simple match/mismatch and substitution-matrix scoring share the same dispatch path. |
| **Affine gaps use Gotoh-style state** | Separate gap-open and gap-extend penalties support realistic scoring. |
| **Named matrix parsing** | Built-in matrix text is parsed into lookup tables. |
| **Codons are packed** | Codon alignment encodes triplets into compact byte tokens. |
| **String wrappers retained** | Strings are converted to typed sequences for compatibility. |

---

## 1. Result Type

### `PairwiseAlignmentResult`

```julia
struct PairwiseAlignmentResult{A <: BioAlphabet} <: AbstractAnalysisResult
    left::BioSequence{A}
    right::BioSequence{A}
    score::Int
    matches::Int
    identity::Float64
    metadata::Dict{Symbol,Any}
end
```

**Description:** Alignment result containing aligned sequences, score, exact match count, identity fraction, and provenance metadata.

---

## 2. Scoring Types

```julia
abstract type AbstractPairwiseScoring end

struct LinearPairwiseScoring <: AbstractPairwiseScoring
    match::Int
    mismatch::Int
end

struct MatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::SubstitutionMatrix
end

struct CodonMatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::CodonSubstitutionMatrix
end
```

**Description:** Internal scoring wrappers used by alignment dispatch.

---

## 3. Substitution Matrices

### `SubstitutionMatrix`

```julia
SubstitutionMatrix(alphabet; match=1, mismatch=-1, default=mismatch, threaded=true)
SubstitutionMatrix(alphabet, scores; default=0, threaded=true)
substitution_matrix(...)
```

**Description:** General matrix container with byte lookup tables and case-folded symbol lookup.

### Named matrices

```julia
named_substitution_matrix(name)
available_named_substitution_matrices()
```

**Description:** Loads built-in substitution matrices such as BLOSUM/PAM-style matrices from `substitution_matrices_data.jl`.

### Codon matrices

```julia
codon_substitution_matrix(alphabet, scores; default=0, threaded=true)
named_codon_substitution_matrix(name)
```

**Description:** Builds or loads 64-codon substitution matrices. Codon alphabets must contain exactly 64 triplets.

---

## 4. Pairwise Alignment

### `pairwise_align`

```julia
pairwise_align(left, right; is_local=false, match=1, mismatch=-1, gap=-1,
               gap_open=nothing, gap_extend=nothing, substitution_matrix=nothing)
```

**Description:** Main pairwise alignment API. Uses global alignment by default and local alignment when `is_local=true`.

**Inputs:** `BioSequence`, byte vectors, strings, `SeqRecordLite`, or `FastqRecord` where supported.

**Gap behavior:**

- `gap` alone gives linear gap scoring;
- `gap_open` or `gap_extend` triggers affine-gap alignment;
- missing affine values default from `gap`.

### Convenience wrappers

```julia
needleman_wunsch(left, right; kwargs...)
smith_waterman(left, right; kwargs...)
local_align(left, right; kwargs...)
```

**Description:** Named wrappers around `pairwise_align`.

---

## 5. Codon Alignment

### `pairwise_align_codons`

```julia
pairwise_align_codons(left, right; is_local=false, kwargs...)
```

**Description:** Aligns DNA sequences in codon triplet space. Sequence lengths must be divisible by 3.

### Codon wrappers

```julia
needleman_wunsch_codons(left, right; kwargs...)
smith_waterman_codons(left, right; kwargs...)
local_align_codons(left, right; kwargs...)
```

**Description:** Named wrappers for codon global/local alignment.

---

## 6. Internal Algorithms

| Helper family | Purpose |
|---|---|
| `_pairwise_align_global` | Linear-gap Needleman-Wunsch. |
| `_pairwise_align_local` | Linear-gap Smith-Waterman. |
| `_pairwise_align_affine_global` | Gotoh affine global alignment. |
| `_pairwise_align_affine_local` | Gotoh affine local alignment. |
| `_pairwise_traceback` | Linear traceback. |
| `_pairwise_traceback_affine` | Affine traceback. |
| `_encode_codon_sequence` | Pack DNA triplets as codon tokens. |
| `_pairwise_align_codon_*` | Codon-space alignment implementations. |

---

## Quick Reference

| API | Purpose |
|---|---|
| `PairwiseAlignmentResult` | Alignment result object. |
| `SubstitutionMatrix` | General substitution matrix. |
| `named_substitution_matrix` | Built-in named matrix. |
| `available_named_substitution_matrices` | List built-in matrices. |
| `codon_substitution_matrix` | Custom codon matrix. |
| `named_codon_substitution_matrix` | Built-in codon matrix. |
| `pairwise_align` | Main DNA/RNA/protein alignment API. |
| `needleman_wunsch` | Global alignment wrapper. |
| `smith_waterman` | Local alignment wrapper. |
| `pairwise_align_codons` | Codon-space alignment API. |

---

## Complete Usage Example

```julia
left = DNASeq("ACGTACGT")
right = DNASeq("ACGTCGT")

global_result = needleman_wunsch(left, right; match=1, mismatch=-1, gap=-1)
local_result = smith_waterman(left, right; match=2, mismatch=-1, gap=-2)

blosum62 = named_substitution_matrix("BLOSUM62")
protein_result = pairwise_align(
    AASeq("MTEYK"),
    AASeq("MTEFK");
    substitution_matrix=blosum62,
    gap=-4)

codon_result = needleman_wunsch_codons(DNASeq("ATGGCC"), DNASeq("ATGACC"))
```
