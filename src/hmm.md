# `hmm.jl` - Hidden Markov Models and Profile HMMs

## Overview

`hmm.jl` implements discrete Hidden Markov Models for biological sequence analysis. It includes Viterbi decoding, forward/backward probabilities, posterior decoding, Baum-Welch training, hard Viterbi training, segmentation, profile HMM construction/scoring, pair-HMM alignment, and lightweight save/load support.

### Purpose

HMMs are a core model family for gene finding, domain detection, sequence segmentation, and profile-based homology search. This module provides those algorithms in pure Julia using log-space probabilities for numerical stability and typed `BioSequence` overloads for BioToolkit integration.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Core probabilities are stored in log space** | Long biological sequences underflow quickly in probability space; log-space dynamic programming is stable. |
| **Discrete byte alphabets** | Emissions are indexed by `UInt8` symbols, matching BioToolkit sequence storage. |
| **Lookup arrays avoid dictionary use in inner loops** | `_emission_lookup` maps byte values to emission columns in O(1). |
| **Training mutates model matrices in place** | `baum_welch!` and `viterbi_train!` preserve the model object while replacing parameters. |
| **Profile HMMs are pragmatic approximations** | The profile builder estimates per-column match/insert/delete parameters from aligned columns with pseudocounts. |
| **I/O avoids a JSON dependency** | `save_hmm` writes a simple key-value text format that `load_hmm` can reconstruct. |

---

## Table of Contents

1. [Core HMM Type](#1-core-hmm-type)
2. [Decoding and Likelihood](#2-decoding-and-likelihood)
3. [Posterior Decoding](#3-posterior-decoding)
4. [Training](#4-training)
5. [Segmentation](#5-segmentation)
6. [Profile HMMs](#6-profile-hmms)
7. [Pair HMMs](#7-pair-hmms)
8. [SummarizedExperiment Integration](#8-summarizedexperiment-integration)
9. [Save and Load](#9-save-and-load)
10. [Quick Reference](#10-quick-reference)

---

## 1. Core HMM Type

### `HMM{T <: AbstractFloat}`

```julia
struct HMM{T<:AbstractFloat}
    states::Vector{String}
    alphabet::Vector{UInt8}
    initial::Vector{T}
    transitions::Matrix{T}
    emissions::Matrix{T}
    _emission_lookup::Vector{Int}
    _unknown_emission::T
end
```

Represents a discrete-emission HMM.

| Field | Description |
|---|---|
| `states` | State names. |
| `alphabet` | Emitted symbols as bytes. |
| `initial` | Initial state log probabilities. |
| `transitions` | State transition log probabilities, indexed `[from, to]`. |
| `emissions` | Emission log probabilities, indexed `[state, alphabet_index]`. |
| `_emission_lookup` | Internal byte-to-column lookup. |
| `_unknown_emission` | Log probability used for symbols not present in the alphabet. |

### Constructor

```julia
HMM(states, alphabet, initial, transitions, emissions; log_space=false)
```

If `log_space=false`, probability-space inputs are converted with `log.(max.(x, eps()))`. If `log_space=true`, the arrays are used as already-log-transformed values.

Validation:

- `length(states) == length(initial)`
- `transitions` is `N x N`
- `emissions` is `N x V`

---

## 2. Decoding and Likelihood

### `viterbi(hmm, sequence) -> (path, log_prob)`

Computes the globally most likely state path.

Supported sequence inputs:

- `AbstractVector{UInt8}`
- `BioSequence`
- `String`

Empty sequences return an empty path and `-Inf`.

### `forward(hmm, sequence) -> log_probability`

Computes total log likelihood of a sequence using the forward algorithm. The implementation uses a rolling two-column buffer for memory efficiency.

### `backward(hmm, sequence) -> Matrix`

Computes the full backward table with dimensions `N x L`. This is used by posterior decoding and Baum-Welch training.

### `hmm_log_likelihood(hmm, sequences) -> Float64`

Sums `forward(hmm, seq)` across a collection of sequences.

---

## 3. Posterior Decoding

### `posterior_state_probabilities(hmm, sequence) -> Matrix{Float64}`

Computes posterior state probabilities:

```text
gamma_t(i) = P(state = i at position t | sequence, hmm)
```

The returned matrix has dimensions `N x L`.

### `posterior_decode(hmm, sequence) -> (assignments, gamma)`

Assigns each position independently to the highest-posterior state. This is softer than Viterbi decoding and does not enforce a globally valid state path.

---

## 4. Training

### `baum_welch!(hmm, sequences; max_iter=100, tol=1e-4, verbose=false)`

Runs Baum-Welch EM training in place.

Returns:

```julia
Vector{Float64} # log-likelihood history
```

Behavior:

- Converts each sequence to byte observations.
- Accumulates expected initial, transition, and emission counts.
- Normalizes rows and writes updated log probabilities back into the HMM.
- Stops when log-likelihood improvement is smaller than `tol`.

### `baum_welch_train(n_states, alphabet, sequences; max_iter=100, tol=1e-4, seed=1, verbose=false)`

Randomly initializes a new HMM and trains it with Baum-Welch.

Returns:

```julia
(hmm, log_likelihoods)
```

### `viterbi_train!(hmm, sequences; max_iter=50, tol=1e-4)`

Runs hard-assignment EM: decode each sequence with Viterbi, then estimate parameters from the decoded paths. It is faster but less statistically complete than Baum-Welch.

---

## 5. Segmentation

### `HMMSegment`

```julia
struct HMMSegment
    start::Int
    stop::Int
    state_index::Int
    state_name::String
    log_probability::Float64
end
```

Represents a contiguous run assigned to one state.

### `segment_sequence(hmm, sequence; method=:viterbi)`

Segments a sequence into state runs.

Methods:

| Method | Assignment source |
|---|---|
| `:viterbi` | Globally optimal state path. |
| `:posterior` | Per-position posterior maximum. |

---

## 6. Profile HMMs

### `ProfileHMMNode`

Stores per-column match emissions, insert emissions, and match/insert/delete transition probabilities.

### `ProfileHMM`

```julia
struct ProfileHMM
    name::String
    length::Int
    alphabet::Vector{UInt8}
    nodes::Vector{ProfileHMMNode}
    null_log_odds::Float64
end
```

Represents a HMMER-style profile over aligned sequence columns.

### `build_profile_hmm(msa; alphabet=nothing, pseudocount=1.0, name="profile")`

Builds a profile HMM from a multiple sequence alignment or vector of aligned strings.

Behavior:

- Requires at least one sequence.
- Requires all sequences to have equal length.
- Defaults to the 20 standard amino-acid alphabet.
- Estimates match emissions per column with pseudocounts.
- Uses gap fraction to set simple match-to-delete probabilities.

### `score_profile_hmm(phmm, sequence; return_per_position=false)`

Scores a sequence against a profile HMM using dynamic programming over match, insert, and delete states. Typed overloads are provided for DNA and amino-acid `BioSequence` values.

---

## 7. Pair HMMs

### `PairHMM`

```julia
struct PairHMM
    match_score::Matrix{Float64}
    gap_open::Float64
    gap_extend::Float64
    alphabet::Vector{UInt8}
end
```

Three-state pair-HMM for pairwise sequence alignment.

### `pair_hmm_align(phmm, seq1, seq2) -> (score, alignment1, alignment2)`

Aligns two typed sequences with affine gap penalties. The DP is equivalent in spirit to Gotoh-style pairwise alignment.

### `pair_hmm_score(phmm, seq1, seq2) -> Float64`

Returns only the alignment score.

---

## 8. SummarizedExperiment Integration

### `baum_welch_from_se(se, n_states; assay_name="counts", kwargs...)`

Treats each assay column as one discretized observation sequence, derives a byte alphabet from rounded assay values, and trains an HMM with `baum_welch_train`.

This is a convenience bridge for matrix-like genomic data, not a replacement for carefully designed observation encoding.

---

## 9. Save and Load

### `save_hmm(path, hmm)`

Writes an HMM in a lightweight text format:

```text
# BioToolkit HMM v1
states=...
alphabet=...
initial=...
trans_1=...
emis_1=...
```

### `load_hmm(path) -> HMM`

Loads a model written by `save_hmm` and reconstructs it with `log_space=true`.

---

## 10. Quick Reference

| Function/Type | Purpose |
|---|---|
| `HMM` | Discrete log-space HMM. |
| `viterbi` | Best state path. |
| `forward` | Sequence log likelihood. |
| `backward` | Backward DP table. |
| `posterior_state_probabilities` | Per-position state posterior matrix. |
| `posterior_decode` | Maximum-posterior assignment. |
| `hmm_log_likelihood` | Sum likelihood across sequences. |
| `baum_welch!` | In-place EM training. |
| `baum_welch_train` | Random initialize and train. |
| `viterbi_train!` | Hard-assignment training. |
| `HMMSegment` | Decoded contiguous state segment. |
| `segment_sequence` | Convert state assignments to runs. |
| `ProfileHMM` | Profile model over aligned columns. |
| `build_profile_hmm` | Build profile HMM from MSA. |
| `score_profile_hmm` | Score sequence against profile HMM. |
| `PairHMM` | Pairwise alignment HMM. |
| `pair_hmm_align` | Pair-HMM alignment and traceback. |
| `pair_hmm_score` | Pair-HMM score only. |
| `baum_welch_from_se` | Train from `SummarizedExperiment`. |
| `save_hmm` | Serialize an HMM. |
| `load_hmm` | Load serialized HMM. |

