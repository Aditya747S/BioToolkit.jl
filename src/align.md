# `align.jl` - Pairwise and Codon-Aware Sequence Alignment

## Overview

`align.jl` implements global, local, affine-gap, matrix-scored, and codon-aware pairwise sequence alignment using well-established algorithms from the bioinformatics literature.

## Algorithms

| Algorithm | Complexity | Reference |
|-----------|------------|-----------|
| Needleman-Wunsch (global) | O(mn) | Needleman & Wunsch, JMB 1970 |
| Smith-Waterman (local) | O(mn) | Smith & Waterman, JMB 1981 |
| Gotoh affine-gap (global) | O(mn) | Gotoh, JMB 1982 |
| Gotoh affine-gap (local) | O(mn) | Gotoh, JMB 1982 |
| Banded Needleman-Wunsch | O(kn) | Ukkonen, 1985 |
| Semi-global alignment | O(mn) | Read alignment variant |
| Overlap alignment | O(mn) | Read assembly variant |
| Codon-aware alignment | O(mn/9) | Frame-preserving |
| Pair-HMM posterior | O(mn) | Durbin et al., 1998 |

## Result Types

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

### `PosteriorAlignmentResult`

```julia
struct PosteriorAlignmentResult{A <: BioAlphabet} <: AbstractAnalysisResult
    posterior_matrix::Array{Float32,3}
    expected_accuracy::Float64
    consensus_alignment::PairwiseAlignmentResult{A}
    log_likelihood::Float64
    state_order::Vector{Symbol}
    metadata::Dict{Symbol,Any}
end
```

**Description:** Posterior-decoded pair-HMM alignment with per-position state probabilities and consensus alignment.

### `GraphAlignmentResult`

```julia
struct GraphAlignmentResult{A <: BioAlphabet} <: AbstractAnalysisResult
    graph_path::Vector{Int}
    graph_sequence::BioSequence{A}
    query_alignment::PairwiseAlignmentResult{A}
    node_scores::Vector{Float64}
    metadata::Dict{Symbol,Any}
end
```

**Description:** Alignment of a query sequence to a sequence graph (DAG).

### `ProfileAlignmentResult`

```julia
struct ProfileAlignmentResult <: AbstractAnalysisResult
    path::Vector{Tuple{Int,Int}}
    score::Float64
    posterior_confidence::Vector{Float64}
    metadata::Dict{Symbol,Any}
end
```

**Description:** Profile-profile alignment path with confidence scores.

## Scoring Models

### Linear Scoring

```julia
struct LinearPairwiseScoring <: AbstractPairwiseScoring
    match::Int
    mismatch::Int
end
```

### Substitution Matrix Scoring

```julia
struct MatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::SubstitutionMatrix
end
```

### Codon Substitution Matrix Scoring

```julia
struct CodonMatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::CodonSubstitutionMatrix
end
```

### Pair-HMM Scoring

```julia
struct PairHMMScoring{S<:AbstractFloat} <: AbstractPairwiseScoring
    match_prob::S
    mismatch_prob::S
    gap_open_prob::S
    gap_extend_prob::S
end
```

**Description:** Probabilistic scoring with match/mismatch emission probabilities and gap state transition probabilities. Uses log-space computation for numerical stability.

### Differentiable Scoring

```julia
struct DifferentiableScoring{F<:AbstractFloat} <: AbstractPairwiseScoring
    weights::Vector{F}
    temperature::F
end
```

**Description:** Smooth scoring for differentiable dynamic programming. As temperature → 0, approaches hard Needleman-Wunsch.

## Substitution Matrices

### `SubstitutionMatrix`

```julia
SubstitutionMatrix(alphabet; match=1, mismatch=-1, default=mismatch, threaded=true)
SubstitutionMatrix(alphabet, scores; default=0, threaded=true)
substitution_matrix(...)
```

General matrix container with byte lookup tables and case-folded symbol lookup. Unknown symbols throw errors rather than silently using a fallback score.

### Named Matrices

```julia
named_substitution_matrix(name)
available_named_substitution_matrices()
```

Built-in matrices (BLOSUM62, BLOSUM80, PAM250) parsed from text with thread-safe caching.

### Codon Matrices

```julia
codon_substitution_matrix(alphabet, scores; default=0, threaded=true)
named_codon_substitution_matrix(name)
```

64×64 codon substitution matrices. Alphabets must contain exactly 64 triplets.

## Main Pairwise Alignment API

### `pairwise_align`

```julia
pairwise_align(left, right; is_local=false, match=1, mismatch=-1, gap=-1,
               gap_open=nothing, gap_extend=nothing, substitution_matrix=nothing)
```

**Description:** Main pairwise alignment API supporting:
- `BioSequence` (DNA, RNA, Protein)
- Byte vectors (`Vector{UInt8}`)
- Strings (converted to DNA)
- `SeqRecordLite` and `FastqRecord`

**Gap Behavior:**
- `gap` alone → linear gap scoring
- `gap_open` or `gap_extend` → affine-gap (Gotoh) alignment
- Missing affine values default from `gap`

**Returns:** `PairwiseAlignmentResult{A}`

### Convenience Wrappers

```julia
needleman_wunsch(left, right; kwargs...)      # global, linear gap
smith_waterman(left, right; kwargs...)        # local, linear gap
local_align(left, right; kwargs...)           # alias for smith_waterman
```

## Advanced Alignment Modes

### Banded Alignment

```julia
banded_needleman_wunsch(seq1, seq2; k=50, match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
```

O(kn) global alignment with band width `k`. Falls back to full DP if |m-n| > k.

### Semi-Global Alignment

```julia
semi_global_align(seq1, seq2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
```

Free end-gaps at both sequence boundaries. Suitable for aligning a read to a reference.

### Overlap Alignment

```julia
overlap_align(seq1, seq2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
```

Free end-gaps at the ends of both sequences. Suitable for read overlap detection in assembly.

### Pair-HMM Posterior Alignment

```julia
pairhmm_align(left, right; scoring=PairHMMScoring())
```

Computes posterior state probabilities via forward-backward algorithm and returns a posterior-decoded consensus alignment with expected accuracy.

### Differentiable Alignment

```julia
differentiable_align(seq1, seq2; scoring=DifferentiableScoring([1.0, -1.0, -1.0]))
soft_alignment_score(left, right, scoring)
```

Smooth objective function approaching hard alignment as temperature → 0. Returns score only (no traceback).

### Unified Dispatcher

```julia
pairalign(mode::Symbol, seq1, seq2; kwargs...)
```

Modes: `:global`, `:local`, `:semi_global`, `:overlap`, `:banded`, `:pair_hmm`, `:differentiable`, `:codon`

## Codon Alignment

### `pairwise_align_codons`

```julia
pairwise_align_codons(left, right; is_local=false, match=1, mismatch=-1,
                      gap=-1, gap_open=nothing, gap_extend=nothing,
                      substitution_matrix=nothing)
```

Aligns DNA sequences in codon triplet space. Sequences must have lengths divisible by 3.

### Codon Wrappers

```julia
needleman_wunsch_codons(left, right; kwargs...)
smith_waterman_codons(left, right; kwargs...)
local_align_codons(left, right; kwargs...)
```

## Coordinate Mapping

```julia
seq2ref(aln, seq_pos)    # query position → reference position
ref2seq(aln, ref_pos)    # reference position → query position
seq2aln(aln, seq_pos)    # query position → alignment column
ref2aln(aln, ref_pos)    # reference position → alignment column
aln2seq(aln, aln_pos)    # alignment column → query position
aln2ref(aln, aln_pos)    # alignment column → reference position
```

All positions are 1-based. Returns 0 if position maps to a gap.

## CIGAR Operations

```julia
cigar(left, right) -> String
cigar(aln::PairwiseAlignmentResult) -> String
parse_cigar(cigar_str) -> Vector{Tuple{Char, Int}}
count_matches(aln)
count_mismatches(aln)
count_insertions(aln)
count_deletions(aln)
count_aligned(aln)
```

SAM-compliant CIGAR string generation and parsing.

## Sequence Graph Alignment

```julia
SequenceGraph(nodes, edges; metadata)
align_to_graph(query, graph; kwargs...)
```

Aligns query to best path through a DAG sequence graph.

## Profile HMM Alignment

```julia
AlignmentProfileHMM(alphabet, match_emissions; insert_emissions, transition_probs, metadata)
align_profiles(profile_a, profile_b; gap=-0.1)
```

Profile-profile alignment using emission dot products.

## Visualization

```julia
visualize_alignment_html(res::PairwiseAlignmentResult; title="...")
visualize_posterior_alignment_html(res::PosteriorAlignmentResult; title="...")
visualize_graph_alignment_html(res::GraphAlignmentResult; title="...")
visualize_sequence_graph_html(graph::SequenceGraph; title="...")
visualize_profile_alignment_html(res::ProfileAlignmentResult; title="...")
to_html(result)  # generic
```

Interactive HTML5/Canvas visualizations.

## Complete Usage Examples

### Basic DNA Alignment

```julia
left = DNASeq("ACGTACGT")
right = DNASeq("ACGTCGT")

# Global alignment (Needleman-Wunsch)
global_result = needleman_wunsch(left, right; match=1, mismatch=-1, gap=-1)

# Local alignment (Smith-Waterman)
local_result = smith_waterman(left, right; match=2, mismatch=-1, gap=-2)

# Affine gap penalties
affine_result = pairwise_align(left, right; match=1, mismatch=-1, gap_open=-5, gap_extend=-1)
```

### Protein Alignment with BLOSUM62

```julia
blosum62 = named_substitution_matrix("BLOSUM62")
protein_result = pairwise_align(
    AASeq("MTEYK"),
    AASeq("MTEFK");
    substitution_matrix=blosum62,
    gap=-4)
```

### Codon Alignment

```julia
# Global codon alignment preserving reading frame
codon_result = needleman_wunsch_codons(DNASeq("ATGGCC"), DNASeq("ATGACC"))

# With custom codon substitution matrix
schneider = named_codon_substitution_matrix("SCHNEIDER")
codon_result2 = pairwise_align_codons(
    DNASeq("ATGGCCATGGCC"),
    DNASeq("ATGACCATGACC");
    substitution_matrix=schneider)
```

### Banded Alignment for Long Sequences

```julia
# Fast alignment when sequences are similar length
banded_result = banded_needleman_wunsch(
    DNASeq("A"^1000),
    DNASeq("A"^995);
    k=50)
```

### Pair-HMM Posterior Decoding

```julia
scoring = PairHMMScoring(match_prob=0.97, mismatch_prob=0.03,
                         gap_open_prob=0.02, gap_extend_prob=0.70)
posterior = pairhmm_align(DNASeq("ACGT"), DNASeq("ACGT"); scoring=scoring)

# Access posterior probabilities
posterior.posterior_matrix[:, :, 1]  # match state probabilities
posterior.expected_accuracy
posterior.consensus_alignment
```

### Coordinate Mapping

```julia
aln = needleman_wunsch(DNASeq("ACGT"), DNASeq("ACGTCGT"))
seq2ref(aln, 3)     # query pos 3 → reference pos
aln2ref(aln, 5)     # alignment col 5 → reference pos
cigar(aln)          # "4M3I"
```

### Semi-Global and Overlap

```julia
# Align read to reference (free end gaps on reference)
semi_global_align(DNASeq("ACGT"), DNASeq("TTACGTCC"))

# Find overlap between reads
overlap_align(DNASeq("ACGTACGT"), DNASeq("TACGTCGT"))
```