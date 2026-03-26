# substitution_matrices_data.jl

## Purpose
This file is a data-only source of substitution matrices. It stores the raw text for many common amino-acid and nucleotide scoring tables so the alignment code can parse and reuse them without hard-coding large matrix literals into algorithmic modules.

## What it contains
- Named amino-acid matrices such as BENNER, BLOSUM, and PAM-style tables.
- Nucleotide matrices such as BLASTN.
- Protein-specific tables such as BLASTP and related scoring variants.
- Extended alphabet tables for ambiguous residues and stop symbols.

## How it is used
Alignment code reads these raw string constants and turns them into numeric scoring matrices when a caller asks for a named substitution scheme. That keeps the data centralized, avoids duplication, and makes it easier to add or audit matrix definitions later.

## Why it matters
Substitution matrices are core inputs for pairwise and local alignment. Keeping them in a dedicated data file separates reference data from scoring algorithms and makes the package easier to maintain.
