# `popgen.jl` - Population Genetics

## Overview

`popgen.jl` provides population/locus containers, allele and genotype frequencies, heterozygosity, HWE tests, F/G statistics, migration estimates, LD, phase inference, genetic distances, PCA/PCoA, effective population size, neutrality/selection statistics, GenePop I/O, Wright-Fisher simulation, GRM/LMM helpers, and coalescent simulation.

### Purpose

This page is a hand-authored reference for `popgen.jl`, grouped around the exported user workflows. Internal helper functions are omitted unless the module exposes them as part of the public API.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Source-matched APIs** | Entries use names implemented in the corresponding `.jl` file. |
| **Workflow sections** | Related functions are documented together so the analysis path is clear. |
| **Concrete result descriptions** | Structs and result types are described by their role in downstream workflows. |
| **Julia-first data flow** | APIs compose through normal Julia arrays, tables, graphs, and BioToolkit objects. |
| **Export-oriented coverage** | The page covers the public functions users are likely to import from BioToolkit. |

---

## 1. Core Types and Frequencies

Population containers hold loci, individuals, and genotypes for classical population genetics.

| API | Description |
|---|---|
| `Locus` | Locus metadata and allele definitions. |
| `PopGenIndividual` | Individual/sample with genotypes and metadata. |
| `Population` | Collection of individuals over loci. |
| `allele_frequencies` | Computes allele frequencies at one locus. |
| `genotype_frequencies` | Computes genotype frequencies at one locus. |
| `heterozygosity_observed` | Observed heterozygosity. |
| `heterozygosity_expected` | Expected heterozygosity. |
| `hardy_weinberg_test` | Approximate HWE test. |
| `hardy_weinberg_exact` | Exact HWE test. |
| `ewens_watterson_test` | Ewens-Watterson neutrality test. |

## 2. Structure, LD, and Distances

These functions quantify differentiation, relatedness, and association between loci/populations.

| API | Description |
|---|---|
| `f_statistics` | Computes F-statistics across populations. |
| `g_statistics` | Computes G-statistics across populations. |
| `migration_rate` | Estimates migration from differentiation. |
| `amova` | Analysis of molecular variance. |
| `linkage_disequilibrium` | LD statistics between two loci. |
| `ld_mapping` | Computes LD across locus windows. |
| `ld_decay` | Summarizes LD decay over physical distance. |
| `infer_phase_em` | Infers two-locus haplotype phase by EM. |
| `genetic_distance` | Computes Nei or related population distances. |
| `population_pca` | PCA over population genotype/frequency features. |
| `population_pcoa` | PCoA from a distance matrix. |
| `mantel_test` | Mantel test between distance matrices. |
| `genetic_relationship_matrix` | Builds GRM from populations. |
| `linear_mixed_model_scan` | Runs an LMM scan using genotype/phenotype/GRM inputs. |
| `inbreeding_coefficient` | Computes individual inbreeding coefficient from GRM. |
| `relatedness` | Computes pairwise relatedness from GRM. |

## 3. Diversity, Selection, and Simulation

Alignment-based methods summarize polymorphism and selection signals.

| API | Description |
|---|---|
| `segregating_sites` | Counts segregating alignment columns. |
| `mismatch_distribution` | Computes pairwise mismatch distribution. |
| `nucleotide_diversity` | Computes pi. |
| `watterson_theta` | Computes Watterson theta. |
| `tajimas_d` | Computes Tajima D. |
| `site_frequency_spectrum` | Builds folded or unfolded SFS. |
| `fu_li_d` | Computes Fu and Li D. |
| `fu_li_f` | Computes Fu and Li F. |
| `ehh` | Extended haplotype homozygosity. |
| `ihs` | Integrated haplotype score. |
| `xp_ehh` | Cross-population EHH. |
| `sweepfinder_clr` | SweepFinder-like composite likelihood ratio scan. |
| `f3_statistic` | Three-population f3 statistic. |
| `f4_statistic` | Four-population f4 statistic. |
| `patterson_d` | Patterson D statistic. |
| `wright_fisher_simulation` | Simulates allele-frequency evolution in one population. |
| `wright_fisher_metapopulation` | Simulates metapopulations with migration. |
| `simulate_coalescent` | Simulates a coalescent genealogy and sequence variation. |

## 4. GenePop and External Wrappers

GenePop utilities support import/export and common external program planning.

| API | Description |
|---|---|
| `GenePopRecord` | Parsed GenePop dataset. |
| `read_genepop_record` | Reads a GenePop file into a record. |
| `read_genepop` | Reads GenePop and converts to populations. |
| `split_in_pops` | Splits a GenePop record by population. |
| `split_in_loci` | Splits a GenePop record by locus. |
| `remove_population!` | Removes a population in place. |
| `remove_locus_by_position!` | Removes a locus by index in place. |
| `remove_locus_by_name!` | Removes a locus by name in place. |
| `write_genepop` | Writes populations as GenePop text. |
| `run_genepop` | Runs/plans GenePop external analysis. |
| `run_fastsimcoal` | Runs/plans fastsimcoal simulations. |
| `run_fdist` | Runs/plans FDIST analysis. |

---

## Complete Usage Example

```julia
using BioToolkit

freq = allele_frequencies(population, 1)
hwe = hardy_weinberg_test(population, 1)
pi = nucleotide_diversity(alignment)
td = tajimas_d(alignment)
```

