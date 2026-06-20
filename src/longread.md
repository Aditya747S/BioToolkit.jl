# `longread.jl` - Long-Read Sequencing

## Overview

`longread.jl` supports Nanopore/PacBio read containers, minimizer sketches, candidate-region search, banded alignment, structural variant evidence, long-read wrappers, phasing, OLC summaries, and chromatin contact extraction.

### Purpose

This page documents the public API in `longread.jl`, grouped by the biological workflow it supports. Each entry describes the role of the function or type in practical BioToolkit analyses.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Domain-specific names** | Public functions match familiar analysis concepts. |
| **Compact result containers** | Typed results preserve the fields needed for downstream inspection. |
| **Composable tables and matrices** | APIs accept common Julia data structures and BioToolkit containers. |
| **Deterministic summaries** | Statistical summaries are reproducible unless stochastic behavior is explicitly requested. |
| **Workflow interoperability** | Outputs can feed plotting, enrichment, clinical, or systems biology modules. |

---

## 1. Types

Containers represent reads, seeds, candidate mapping intervals, alignments, and structural variants.

| API | Description |
|---|---|
| `NanoporeRead` | Nanopore read sequence with quality and metadata. |
| `PacBioRead` | PacBio read sequence with quality and metadata. |
| `StructuralVariant` | Structural variant call with type, coordinates, evidence, and score. |
| `MinimizerSeed` | Minimizer hash and position used for seeding. |
| `CandidateRegion` | Candidate mapping interval from minimizer hits. |
| `BandedAlignmentResult` | Banded semiglobal alignment score, path, and coordinates. |

## 2. Sketching and Alignment

Minimizer and banded-alignment helpers provide a pure-Julia mapping skeleton.

| API | Description |
|---|---|
| `canonical_kmer_hash` | Computes strand-canonical k-mer hashes. |
| `minimizer_sketch` | Builds minimizer sketches for long reads or references. |
| `build_minimizer_index` | Indexes reference minimizers. |
| `find_candidate_regions` | Finds likely alignment regions from minimizer hits. |
| `banded_semiglobal_alignment` | Aligns a read to a candidate region within a diagonal band. |
| `minimap2_align` | Wrapper/payload helper for minimap2-style alignment. |
| `graphaligner_align` | Wrapper/payload helper for graphaligner-style alignment. |

## 3. Structural Variation and Long-Read Workflows

Evidence functions summarize split reads, discordant pairs, phasing, contacts, and read metrics.

| API | Description |
|---|---|
| `split_read_evidence` | Extracts split-read evidence for breakpoints. |
| `discordant_pair_evidence` | Summarizes discordant read-pair evidence. |
| `cluster_structural_variants` | Clusters breakpoint evidence into SV candidates. |
| `call_structural_variants` | Runs the built-in SV calling workflow. |
| `sniffles_call` | Wrapper/payload helper for Sniffles-like SV calling. |
| `svim_call` | Wrapper/payload helper for SVIM-like SV calling. |
| `isoseq_consensus` | Builds Iso-Seq transcript consensus summaries. |
| `whatshap_phase` | Wrapper/payload helper for WhatsHap-style phasing. |
| `phase_reads_by_alleles` | Assigns reads to haplotypes using allele observations. |
| `overlap_layout_consensus` | Builds overlap-layout-consensus summaries. |
| `read_n50` | Computes read N50. |
| `overlap_graph_table` | Creates read overlap graph tables. |
| `porec_contact_table` | Extracts Pore-C contact records. |
| `omnic_contact_matrix` | Builds Omni-C contact matrices. |

---

## Complete Usage Example

```julia
using BioToolkit

sketch = minimizer_sketch(read_sequence; k=15, w=10)
index = build_minimizer_index(reference_sequences)
regions = find_candidate_regions(sketch, index)
aln = banded_semiglobal_alignment(read_sequence, reference_window)
```

