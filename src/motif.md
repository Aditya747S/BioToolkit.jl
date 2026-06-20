# `motif.jl` - Sequence Motif Discovery, Scanning, and Motif I/O

## Overview

`motif.jl` implements BioToolkit's motif model stack: count matrices, frequency matrices, position weight matrices, scanning, de novo seed-based discovery, sequence-logo SVG rendering, and parsers for common motif formats. The API is typed around `BioSequence{A}` and `BioAlphabet`, so nucleotide and amino-acid workflows can share the same machinery while still preserving alphabet dispatch.

### Purpose

Motif analysis needs three things to be reliable: a clear symbol alphabet, a reproducible scoring model, and interoperable import/export of motif databases. This module provides those pieces without requiring external motif libraries.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Motif matrices carry alphabet order** | Rows are interpreted by `alphabet`, avoiding hard-coded assumptions about DNA-only matrices. |
| **PWM values are probabilities** | `MotifPWM` stores normalized position probabilities; scanning converts them to log2 odds against a background model. |
| **Typed hits preserve window alphabet** | `MotifHit{A}` and `MotifSite{A}` keep the matched sequence as a typed `BioSequence{A}`. |
| **Backgrounds are normalized internally** | `nothing`, vectors, and dictionaries are accepted and converted to probability vectors. |
| **Reverse-strand scanning is nucleotide-only** | `motif_scan_both_strands` is restricted to DNA/RNA alphabets, where reverse complements are meaningful. |
| **Format parsers produce `MotifProfile`** | MEME, AlignACE, and JASPAR imports retain metadata, counts, PWM, and occurrence information where available. |

---

## Table of Contents

1. [Matrix Types](#1-matrix-types)
2. [Hit and Discovery Types](#2-hit-and-discovery-types)
3. [Counting and Matrix Construction](#3-counting-and-matrix-construction)
4. [Information Metrics](#4-information-metrics)
5. [Scanning](#5-scanning)
6. [De Novo Discovery](#6-de-novo-discovery)
7. [Motif Profiles and I/O](#7-motif-profiles-and-io)
8. [Sequence Logos](#8-sequence-logos)
9. [Quick Reference](#9-quick-reference)

---

## 1. Matrix Types

### `MotifCounts{A}`

```julia
struct MotifCounts{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    counts::Matrix{Int}
end
```

Stores raw per-position symbol counts. Rows correspond to `alphabet`; columns correspond to motif positions.

### `MotifFrequencyMatrix{A}`

```julia
struct MotifFrequencyMatrix{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    values::Matrix{Float64}
end
```

Stores per-column probabilities. It is usually produced from `MotifCounts` with optional pseudocounts.

### `MotifPWM{A}`

```julia
struct MotifPWM{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    values::Matrix{Float64}
end
```

Stores position probabilities used for log2-odds scoring during scans. Despite the name, values are not precomputed log odds; background correction happens at scan time.

Each matrix type has a constructor accepting `Char` or `UInt8` alphabets and numeric matrices. The alphabet string is used to infer a BioToolkit alphabet type.

---

## 2. Hit and Discovery Types

### `MotifHit{A}`

```julia
struct MotifHit{A <: BioAlphabet}
    start::Int
    strand::Int8
    score::Float64
    window::BioSequence{A}
end
```

Represents one PWM scan hit.

| Field | Description |
|---|---|
| `start` | One-based start coordinate in the scanned sequence. |
| `strand` | `1` for forward, `-1` for reverse complement. |
| `score` | Log2-odds motif score. |
| `window` | Matched sequence window. |

### `MotifSite{A}`

Stores one discovered motif site, including the source sequence index, start coordinate, strand, mismatch count, and oriented window.

### `MotifDiscoveryResult{A}`

```julia
struct MotifDiscoveryResult{A <: BioAlphabet} <: AbstractAnalysisResult
    seed::BioSequence{A}
    alphabet::Vector{UInt8}
    counts::MotifCounts{A}
    pwm::MotifPWM{A}
    sites::Vector{MotifSite{A}}
    support::Int
    information_content::Float64
end
```

Result object returned by `discover_motifs`. It integrates with BioToolkit's analysis-result summary machinery.

---

## 3. Counting and Matrix Construction

### `motif_counts(sequences; alphabet="ACGT")`

Counts aligned motif instances.

Supported inputs:

- `AbstractVector{BioSequence{A}}`
- `AbstractVector{<:AbstractString}` (converted to DNA sequences)
- `AbstractVector{<:SeqRecordLite}`

Validation:

- At least one sequence is required.
- All sequences must have the same length.
- Every observed symbol must be in the requested alphabet.

Example:

```julia
sites = DNASeq.(["ACGT", "ACGA", "ACGG"])
counts = motif_counts(sites; alphabet="ACGT")
```

### `motif_frequency_matrix(profile; pseudocount=0.0)`

Converts a count matrix to per-column probabilities. `pseudocount` must be nonnegative.

### `motif_pwm(profile; pseudocount=0.0, background=nothing)`

Builds a `MotifPWM` from counts or from sequences/records. When a background is supplied, the pseudocount is distributed in proportion to background frequencies.

Supported background forms:

| Background | Meaning |
|---|---|
| `nothing` | Uniform over the motif alphabet. |
| `AbstractVector` | Values are normalized to sum to one. |
| `AbstractDict` | Keys may be `Char` or `UInt8`; values are normalized. |

### `motif_consensus(profile; threshold=0.5)`

Returns a consensus string from count columns. A column becomes `N` when the best symbol is below the threshold, tied, or unsupported.

---

## 4. Information Metrics

### `motif_entropy(profile::MotifFrequencyMatrix)`

Computes average Shannon entropy per motif column.

### `motif_relative_entropy(profile::MotifFrequencyMatrix; background=nothing)`

Computes average relative entropy against a background distribution.

### `motif_information_content(profile::MotifPWM; background=nothing)`

Computes average motif information content from PWM probabilities.

All background-aware methods require positive background entries for any symbol with nonzero motif probability.

---

## 5. Scanning

### `motif_scan(sequence::BioSequence{A}, pwm::MotifPWM{A}; threshold=0.0, background=nothing)`

Scans a typed sequence on the forward strand.

Algorithm:

1. Build a byte-to-row lookup from `pwm.alphabet`.
2. Slide a window of motif width along the sequence.
3. For each valid window, sum `log2(probability / background_probability)`.
4. Return windows whose score is at least `threshold`.

Invalid symbols in a window cause that window to be skipped.

### `motif_scan_both_strands(sequence::BioSequence{A}, pwm::MotifPWM{A}; threshold=0.0, background=nothing)`

Scans both forward and reverse-complement windows for DNA/RNA sequences. Returned hits are sorted by `(start, -strand, -score)`.

Example:

```julia
counts = motif_counts(DNASeq.(["ACGT", "ACGA", "ACGG"]))
pwm = motif_pwm(counts; pseudocount=0.5)
hits = motif_scan_both_strands(DNASeq("TTACGTAA"), pwm; threshold=0.0)
```

---

## 6. De Novo Discovery

### `discover_motifs(sequences; kwargs...)`

Performs seed-and-extend motif discovery.

Key options:

| Option | Default | Description |
|---|---:|---|
| `alphabet` | `"ACGT"` for DNA, amino-acid alphabet otherwise | Symbol order for count/PWM construction. |
| `k` | `6` | Motif seed length. |
| `top_n` | `3` | Maximum number of motifs returned. |
| `min_support` | `2` | Minimum number of sequences/sites supporting a motif. |
| `max_mismatches` | `1` | Maximum mismatches allowed when assigning a site to a seed. |
| `pseudocount` | `0.5` | PWM smoothing value. |
| `background` | `nothing` | Background distribution. |
| `reverse_complements` | `true` | Canonicalize DNA seeds with reverse complements. |

The method ranks frequent k-mers, finds the best compatible site per sequence, builds counts and a PWM from selected sites, and returns `MotifDiscoveryResult` values.

---

## 7. Motif Profiles and I/O

### `MotifOccurrence{A}`

Represents one motif occurrence from imported motif files.

Fields:

- `sequence_id`
- `start`
- `strand`
- `score`
- `pvalue`
- `site`

### `MotifProfile`

```julia
struct MotifProfile
    name::String
    alphabet::Vector{Char}
    counts::MotifCounts
    pwm::MotifPWM
    occurrences::Vector{MotifOccurrence}
    metadata::Dict{String,String}
end
```

Format parser output that keeps the named motif, count matrix, PWM, occurrences, and source metadata together.

### `read_meme(path_or_io; alphabet="ACGT")`

Reads MEME letter-probability matrices and optional occurrence-like lines. Returns `Vector{MotifProfile}`.

### `read_alignace(path_or_io; alphabet="ACGT")`

Reads AlignACE-style motif blocks, collects sites, builds counts and smoothed PWM values, and returns motif profiles.

### `read_jaspar(path_or_io; alphabet="ACGT")`

Reads JASPAR-style matrices into motif profiles.

The parsers are provenance-aware when an active context exists and store source/hash metadata where available.

---

## 8. Sequence Logos

### `sequence_logo_svg(profile; width=720, height=220, background=nothing, title=nothing)`

Returns an SVG string showing per-position motif information. The implementation stacks colored rectangles and letters by probability contribution.

Color defaults:

| Letter | Color |
|---|---|
| `A` | Green |
| `C` | Blue |
| `G` | Orange |
| `T`/`U` | Red |
| Other | Gray |

### `motif_logo_svg(profile; kwargs...)`

Alias for `sequence_logo_svg`.

### `sequence_logo(profile; kwargs...)`

Alias for `sequence_logo_svg`.

---

## 9. Quick Reference

| Function/Type | Purpose |
|---|---|
| `MotifCounts` | Raw count matrix. |
| `MotifFrequencyMatrix` | Probability matrix. |
| `MotifPWM` | Scanning PWM probability matrix. |
| `MotifHit` | PWM scan hit. |
| `MotifSite` | De novo motif site. |
| `MotifDiscoveryResult` | De novo discovery result. |
| `MotifProfile` | Imported motif profile. |
| `motif_counts` | Count aligned motif sites. |
| `motif_frequency_matrix` | Convert counts to frequencies. |
| `motif_pwm` | Build PWM probabilities. |
| `motif_entropy` | Average column entropy. |
| `motif_relative_entropy` | Average relative entropy. |
| `motif_information_content` | Average information content. |
| `motif_scan` | Forward-strand PWM scan. |
| `motif_scan_both_strands` | DNA/RNA forward and reverse-complement scan. |
| `discover_motifs` | Seed-based motif discovery. |
| `motif_consensus` | Consensus sequence from counts. |
| `sequence_logo_svg` | SVG motif logo. |
| `read_meme` | Parse MEME motifs. |
| `read_alignace` | Parse AlignACE motifs. |
| `read_jaspar` | Parse JASPAR motifs. |

