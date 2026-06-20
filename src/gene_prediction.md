# `gene_prediction.jl` - Ab Initio Gene Prediction

## Overview

`gene_prediction.jl` implements lightweight gene-structure analysis for DNA sequences: HMM coding-region detection, start/stop codon discovery, Kozak scoring, splice-site detection, coding-potential metrics, GFF3 export, and prediction summaries.

### Purpose

Use this module when a sequence needs first-pass gene prediction without an external annotation pipeline. It favors inspectable heuristics and simple data structures over opaque model files.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Domain-specific containers** | Results use explicit structs where the workflow has stable fields or needs provenance metadata. |
| **Workflow APIs** | Functions expose complete analysis steps, not only internal numeric kernels. |
| **Table-compatible outputs** | Outputs are designed to work with Julia arrays, dictionaries, named tuples, and DataFrames. |
| **Pure-Julia core** | External tools are optional wrappers; core summaries remain usable inside Julia. |
| **Provenance hooks** | Many public functions accept provenance context keywords or return result containers with provenance records. |

---

## 1. Core Types

Prediction outputs are represented with interval-level containers that can be exported or attached to annotated records.

| API | Description |
|---|---|
| `GeneInterval` | Minimal `(start, stop)` interval returned by the HMM coding-region detector. |
| `ExonInterval` | Predicted exon span with strand, frame phase, end phase, and score. |
| `IntronInterval` | Predicted intron span with donor and acceptor splice-site scores. |
| `SpliceSite` | Detected donor or acceptor motif with position, strand, score, and matched consensus. |
| `GenePrediction` | Complete gene model with gene id, chromosome, strand, exons, introns, CDS/protein length, score, and source. |

## 2. Signals and Coding Scores

These functions identify sequence signals that feed gene models or quality-control reports.

| API | Description |
|---|---|
| `find_start_codons` | Finds ATG start codons on both strands and annotates them with Kozak scores. |
| `find_stop_codons` | Finds TAA/TAG/TGA stop codons across all reading frames on both strands. |
| `score_kozak_context` | Scores the nucleotide context around an ATG using a Kozak-style position-weight table. |
| `detect_splice_sites` | Detects donor and acceptor splice motifs with score filtering. |
| `calculate_testcode` | Computes Fickett/TestCode-style coding potential in sliding windows. |
| `codon_bias_index` | Scores codon-bias support against an optional reference codon table. |

## 3. Prediction and Reporting

High-level functions assemble, filter, summarize, and export predictions.

| API | Description |
|---|---|
| `predict_genes_hmm` | Uses a two-state noncoding/coding HMM and Viterbi decoding to return `GeneInterval`s. |
| `predict_gene_structure` | Combines ORF, splice-site, and coding-potential evidence into `GenePrediction`s. |
| `predict_utr_regions` | Adds approximate UTR regions around predicted coding models. |
| `gene_density` | Computes gene density over a total sequence/genome length. |
| `gff3_export` | Exports predictions in GFF3-style feature text. |
| `cds_statistics` | Summarizes CDS and protein lengths across predictions. |
| `filter_gene_predictions` | Filters predictions by score, length, exon count, and related criteria. |
| `gene_prediction_summary` | Builds a compact aggregate summary of prediction counts and lengths. |
| `annotate_orf_features` | Adds predicted ORF/gene features to annotated sequence records. |

---

## Quick Reference

| Area | Main APIs |
|---|---|
| Core Types | `GeneInterval`, `ExonInterval`, `IntronInterval`, `SpliceSite`, `GenePrediction` |
| Signals and Coding Scores | `find_start_codons`, `find_stop_codons`, `score_kozak_context`, `detect_splice_sites`, `calculate_testcode`, `codon_bias_index` |
| Prediction and Reporting | `predict_genes_hmm`, `predict_gene_structure`, `predict_utr_regions`, `gene_density`, `gff3_export`, `cds_statistics`, ... |

---

## Complete Usage Example

```julia
using BioToolkit

seq = DNASeq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
starts = find_start_codons(seq)
stops = find_stop_codons(seq)
splice = detect_splice_sites(seq)
intervals = predict_genes_hmm(seq)
summary = gene_density(intervals, length(seq); window_bp=1000)
```

