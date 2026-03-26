# motif.jl

## Purpose
This file implements motif discovery, motif counting, and motif-scoring utilities for DNA sequence analysis. It provides both the intermediate representations and the scanning routines needed for discovery-style and annotation-style motif work.

## Main types
- `MotifCounts` stores per-position base counts.
- `MotifFrequencyMatrix` stores normalized position frequencies.
- `MotifPWM` stores the final scoring matrix used for scanning.
- `MotifHit` stores a motif match on a sequence window.
- `MotifSite` stores motif site information across many sequences.
- `MotifDiscoveryResult` bundles a discovered seed, counts, PWM, sites, support, and information content.

## Public functions
- Counting and matrix building: `motif_counts` and `motif_frequency_matrix`.
- Scoring summaries: `motif_entropy`, `motif_relative_entropy`, and `motif_information_content`.
- Scanning and discovery: `motif_scan` and the discovery helpers built around motif seeds.
- Convenience helpers: `DataFrame(::Vector{MotifHit})` and `DataFrame(::Vector{MotifSite})` conversions.

## How it is used
The usual workflow is to count aligned candidate sites with `motif_counts`, convert them into a `MotifPWM` through `motif_frequency_matrix`, then score a larger sequence with `motif_scan`.

The discovery helpers support seed-based enumeration and reverse-complement-aware site selection, which makes the module useful for finding enriched k-mers in regulatory sequence sets.

## Implementation notes
- The module uses compact byte lookup tables so common nucleotide characters can be handled efficiently.
- Pseudocounts are supported in the frequency-matrix builder to avoid zero-probability columns.
- The `DataFrame` conversions make it easy to turn motif hits or sites into tabular outputs for downstream filtering or plotting.

## Why it matters
Motif analysis connects raw sequences to regulatory interpretation. This file gives BioToolkit the basic building blocks needed to discover motifs, score them, and present them in a table-friendly way.
