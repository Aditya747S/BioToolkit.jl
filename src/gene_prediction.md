# gene_prediction.jl

## Purpose
This file implements a compact hidden Markov model gene predictor. It is a lightweight ab initio gene-finding utility that scans a nucleotide sequence, classifies bases into coding and non-coding states, and returns the likely gene intervals.

## Main function
- predict_genes_hmm(sequence; p_coding=0.01): run a two-state HMM over a DNA sequence and return an array of (start_index, end_index) tuples for predicted genes.

## What the function does
The function builds a small HMM with two states: non-coding background and coding sequence. It uses a simple nucleotide alphabet, fixed initial probabilities, a transition matrix that favors long coding segments, and emission probabilities that prefer GC-rich coding regions. The sequence is decoded with Viterbi, then contiguous coding runs are extracted as candidate genes.

## Internal workflow
1. Encode the DNA sequence as bytes.
2. Define the HMM states, alphabet, priors, transitions, and emissions.
3. Run Viterbi to obtain the most likely state path.
4. Convert stretches of coding-state positions into genomic intervals.
5. Filter out short, low-confidence segments below the minimum length threshold.

## Output
The function returns a vector of tuples, each tuple containing the start and end coordinates of a predicted gene. The coordinates are one-based sequence positions, making them easy to feed into downstream annotation, visualization, or export steps.

## How it is used
This module is most useful when a user has raw sequence and wants a quick, model-based pass before deeper annotation. It can be used to suggest candidate ORFs or to provide rough coding-region boundaries that can later be refined by other tools.

## Why this file matters
This file demonstrates how BioToolkit uses its generic HMM machinery for a concrete biological task. It is intentionally small, but it gives the package a simple gene-prediction capability that does not require an external pipeline.
# Gene Prediction API (gene_prediction.jl)

The native applied gene parsing module ties directly into our robust CPU Profile HMM `viterbi` framework to solve genomic overlap and naive false-positive constraints common to greedy `find_orfs` heuristic sweeps.

### Predict Genes (HMM Abstract)
The `predict_genes_hmm(sequence; p_coding=0.01)` constructs a mathematical 2-State (Coding vs NonCoding) boundary model isolating elevated genomic GC content spikes bounded strictly by transitional density probabilities.

```julia
using BioToolkit

# Pass raw raw sequence fragments (it inherently blocks meaningless chunks ≤ 30bp)
genes = predict_genes_hmm("ATATAAATTTTAAATATATATAAGCGCGCGCGCGCGCGCGCGCGCGCATA")

# Emits tuple boundaries [start, stop] representing Exon blocks
println(genes) # e.g. [(24, 47)]
```

### Profile Structural Architecture
Internally, the pipeline allocates a localized statistical machine avoiding hardcoded substring targets.
- **State 1 (NonCoding)**: GC distribution flat (background), representing noise.
- **State 2 (Coding)**: GC distribution artificially elevated mathematically (simulating generalized Exon bounds).
- **Viterbi Trace**: It executes a full Log-Space trace and isolates State 2 continuity. The predictive nature mathematically overrides internal random noise (introns/defects) due to the Log-Sum sequence dependence.
