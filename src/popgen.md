# popgen.jl

## Purpose
This file provides the population-genetics data model and the core frequency/statistics helpers. It is the lightweight base layer for locus-level summaries across individuals and populations.

## Main types
- `Locus{T}` stores the alleles observed at one locus.
- `PopGenIndividual{T}` stores an individual ID and its loci.
- `Population{T}` stores a collection of individuals and the population name.

## Public functions
- `allele_frequencies(pop, locus_idx)` computes allele frequencies for a locus.
- `genotype_frequencies(pop, locus_idx)` computes genotype frequencies using unordered allele sets.
- `heterozygosity_observed(pop, locus_idx)` computes the observed heterozygote proportion.
- `heterozygosity_expected(pop, locus_idx)` computes expected heterozygosity from allele frequencies.
- `hardy_weinberg_test(pop, locus_idx)` computes a fast Hardy-Weinberg equilibrium p-value approximation.
- `hardy_weinberg_exact(pop, locus_idx)` computes an exact Hardy-Weinberg test for bi-allelic loci.
- `f_statistics(populations, locus_idx)` computes F-statistics across subpopulations.

## How it is used
Users build a `Population` from `PopGenIndividual` values, then call the frequency and heterozygosity helpers to summarize variation at a locus. Those summaries feed downstream equilibrium tests and differentiation metrics.

The module is deliberately small and composable: it does not try to cover every population-genetics statistic, but it provides the core data structures that more specialized routines can build on.

## Implementation notes
- The code is written around arbitrary allele element types through parametric structs.
- `hardy_weinberg_test` uses a simplified chi-square approximation for a fast path.
- `hardy_weinberg_exact` handles the exact bi-allelic case explicitly.
- The file-level overview in the source advertises broader population-genetics capabilities elsewhere in the package, but the code in this module itself is centered on the data containers and locus summaries listed above.

## Why it matters
Population genetics needs a clean representation of loci, individuals, and populations before higher-level summaries can be computed. This file supplies that foundation and keeps the basic locus-level calculations in one place.
