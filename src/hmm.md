# hmm.jl

## Purpose
This file implements a general-purpose hidden Markov model engine in log-space. It is designed to be reusable by other modules that need numerically stable sequence inference, such as gene prediction or other probabilistic annotators.

## Main type
- `HMM{T}` stores the state list, alphabet, initial distribution, transition matrix, emission matrix, a byte-to-alphabet lookup table, and a fallback probability for unknown symbols.

## Constructors
- `HMM{T}(states, alphabet, initial, transitions, emissions; unknown_prob)` builds a typed model and validates matrix dimensions.
- `HMM(states, alphabet, initial, transitions, emissions; log_space=false)` is a convenience constructor that converts ordinary probabilities into log-space.

## Public functions
- `viterbi(hmm, sequence)` returns the most likely hidden-state path and its log-probability.
- `forward(hmm, sequence)` returns the total log-probability across all valid paths.
- `backward(hmm, sequence)` returns the backward-probability matrix.

## Internal helpers
- `_get_emission_logprob` resolves the emission log-probability for a state and byte.
- `logaddexp` performs numerically stable log-space addition.

## How it is used
Users create an `HMM` from a state set, an alphabet, and probability tables, then call `viterbi`, `forward`, or `backward` on either a byte sequence or a string. The string overloads convert through `codeunits`, so the model works directly on DNA, protein, or other ASCII-encoded sequence data.

## Implementation notes
- All probabilities are stored in log-space to reduce underflow risk.
- The emission lookup table makes byte-to-column mapping constant time.
- Unknown symbols fall back to the configured unknown emission probability.

## Why it matters
Probabilistic sequence models are a shared need across several BioToolkit modules. This file keeps the HMM implementation centralized, stable, and reusable instead of duplicating inference code in each analysis layer.
