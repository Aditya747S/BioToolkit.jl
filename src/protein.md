# protein.jl

## Purpose
This file implements fast protein property calculations in a table-driven style. It mirrors the kinds of outputs users expect from ExPASy ProtParam without adding heavy dependencies or runtime lookup overhead.

## Main data tables
- `_AA_MASS_MONO` and `_AA_MASS_AVG` store residue masses for monoisotopic and average molecular weight calculations.
- `_HYDROPATHICITY` stores Kyte-Doolittle hydropathy scores.
- `_AA_CODE` and `_AA_CHARS` provide compact amino-acid encoding for fast indexing.
- `_INSTABILITY_DIWV` stores the Guruprasad dipeptide instability weights used by the instability index calculation.
- `_PK_NTERM`, `_PK_CTERM`, `_PK_D`, `_PK_E`, `_PK_C`, `_PK_Y`, `_PK_H`, `_PK_K`, and `_PK_R` hold the pK constants used by the isoelectric-point estimator.

## Public functions
- `protein_mass(sequence; type="monoisotopic")` computes molecular weight in Daltons.
- `extinction_coefficient(sequence)` computes the 280 nm extinction coefficient using the Pace formula.
- `instability_index(sequence)` computes the Guruprasad instability index.
- `isoelectric_point(sequence; precision=0.01)` estimates pI using a bisection-style IPC algorithm.
- `gravy(sequence)` computes the grand average of hydropathicity.
- `aliphatic_index(sequence)` computes the aliphatic index from Ala, Val, Ile, and Leu content.
- `protparam(sequence)` returns a bundled summary of the main physicochemical properties.

## How it is used
The basic functions work on any amino-acid string and are intended for direct use in annotation or reporting pipelines. For example, `protein_mass` is useful for mass validation, `gravy` and `aliphatic_index` help describe hydrophobicity, and `instability_index` provides a coarse in vitro stability signal.

`protparam` is the best entry point when a caller wants a single summary object. It performs one pass over the sequence, then returns a `NamedTuple` containing:
- `length`
- `molecular_weight_mono`
- `molecular_weight_avg`
- `negative_residues`
- `positive_residues`
- `extinction_coefficient`
- `instability_index`
- `aliphatic_index`
- `gravy`
- `isoelectric_point`

## Implementation notes
- Lookups are intentionally done through 256-element arrays so ASCII residue access stays constant time.
- Functions throw `ArgumentError` when they encounter unknown residues or invalid input lengths.
- The module treats protein properties as sequence-level summaries, not as structural predictions.

## Threading notes
- This module is intentionally left serial: the main operations are already O(n) with very small per-residue work, so `Threads.@threads` would add overhead without a meaningful speedup.
- If protein summaries are needed in bulk, parallelize at the caller level over many sequences instead of inside these helpers.

## Why it matters
Protein descriptors are used throughout annotation, proteomics, and downstream interpretation. This file gives BioToolkit a compact, predictable implementation of the standard sequence statistics people reach for most often.
