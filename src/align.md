# align.jl

## Purpose
This file implements the pairwise alignment engine used throughout BioToolkit. It does two jobs at once: it defines the data structures that describe an alignment result, and it implements the scoring and traceback logic that produces those results for nucleotide, protein, and codon-level sequences.

## Main structs
- PairwiseAlignmentResult: the final alignment object, including the aligned left and right strings, the total score, the number of matches, and the identity fraction.
- AbstractPairwiseScoring: the abstract parent type for all scoring models.
- LinearPairwiseScoring: a simple match/mismatch scoring model.
- SubstitutionMatrix: a general substitution matrix container with alphabet lookup tables and a default score.
- MatrixPairwiseScoring: a scoring wrapper around a substitution matrix.
- CodonSubstitutionMatrix: a substitution matrix specialized for codon tokens.
- CodonMatrixPairwiseScoring: the codon-scoring wrapper used by codon alignment routines.

## Important constants and embedded data
- _STANDARD_SUBSTITUTION_MATRIX_TEXT stores built-in amino-acid matrices as raw text, including BLOSUM62, BLOSUM80, and PAM250.
- _STANDARD_SUBSTITUTION_MATRIX_CACHE memoizes parsed matrices so repeated requests do not rebuild the same object.
- _CODON_DECODE_TABLE and the codon encoding helpers support translation-aware alignment.
- Internal traceback state constants track match, gap-left, gap-right, and no-trace states during dynamic programming.

## Public functions
- SubstitutionMatrix(alphabet; match, mismatch, default): builds a simple square substitution matrix from a character alphabet.
- SubstitutionMatrix(alphabet, scores; default): wraps an explicit score matrix.
- named_substitution_matrix(name): loads one of the named built-in matrices.
- named_codon_substitution_matrix(name): loads a codon-level substitution matrix by name.
- codon_substitution_matrix(alphabet, scores; default): constructs a codon-aware scoring matrix from explicit values.
- pairwise_align(...): performs pairwise alignment using the selected scoring scheme and gap model.
- pairwise_align_codons(...): performs codon-aware pairwise alignment.

## What the alignment routines do
The alignment pipeline uses dynamic programming to compare two biological sequences and then reconstructs the best-scoring path. Depending on the function variant, it can run with a simple linear gap model, an affine gap model, or a codon-aware model. The file also includes private helper routines for the traceback phase, because traceback logic differs between global, local, and affine-gap alignments.

## Typical usage
1. Choose or build a scoring model. For example, use named_substitution_matrix("BLOSUM62") for protein alignment or build a simple LinearPairwiseScoring for quick DNA comparisons.
2. Call pairwise_align(left, right; scoring=..., gap=...) or pairwise_align_codons(...) depending on the biological representation.
3. Inspect the returned PairwiseAlignmentResult to see the aligned strings, total score, and identity.
4. Use the result object as input for downstream reporting, visualization, or sequence filtering.

## Internal helpers
The file contains a substantial amount of private machinery, including matrix parsing, codon token packing and decoding, and separate global/local/affine alignment kernels for standard and codon-based sequences. Those helpers are not part of the normal user-facing API, but they are what make the public alignment functions fast and flexible.

## Threading notes
- Substitution-matrix parsing and construction now default to threaded row/entry initialization where the work is independent.
- The core pairwise dynamic-programming kernels are still recurrence-bound, so the main alignment pass remains serial until the recurrence itself is rewritten in a wavefront-parallel form.

## Why this file matters
This is the alignment backbone for the package. Other modules do not need to reinvent scoring or traceback; they can rely on the result type and the alignment routines here.
