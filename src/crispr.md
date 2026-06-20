# `crispr.jl` - CRISPR Guide Design

## Overview

`crispr.jl` designs and scores guide RNAs for Cas systems, base editors, prime editing, HDR templates, genome-wide off-target summaries, and CRISPR screen analysis.

### Purpose

This page documents the user-facing types and workflows implemented in `crispr.jl`. The emphasis is on public APIs exported through BioToolkit and the biological analysis step each one supports.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Biological workflow names** | APIs are named after common analysis tasks so code reads like a methods section. |
| **Typed summaries where useful** | Reusable outputs have named structs instead of anonymous dictionaries. |
| **Matrix/table interoperability** | Functions accept standard Julia matrices, vectors, and table-like records for easy integration. |
| **Deterministic defaults** | Randomized or approximate methods expose seeds, thresholds, or iteration controls where relevant. |
| **Downstream-ready results** | Outputs can feed plotting, enrichment, annotation, or reporting modules. |

---

## 1. Types and Systems

System and result types describe guides, PAMs, off-targets, editing windows, and supported nucleases/editors.

| API | Description |
|---|---|
| `GuideRNA` | Guide sequence, target coordinates, strand, PAM, and scores. |
| `CRISPRSystem` | Cas/editor metadata including PAM rules and guide length. |
| `OffTarget` | Candidate off-target with mismatches, position, and score. |
| `EditingWindow` | Editable base window for base-editor or prime-editor design. |
| `SpCas9` | Canonical NGG SpCas9 system. |
| `SaCas9` | SaCas9 system preset. |
| `Cas12a` | Cas12a/Cpf1-style preset. |
| `CasX` | CasX-style preset. |
| `SpRY` | Relaxed-PAM SpRY preset. |
| `BE3` | C-to-T base-editor preset. |
| `ABE8e` | A-to-G base-editor preset. |
| `CBE4max` | C-to-T base-editor preset. |

## 2. Guide Discovery and Scoring

These functions scan sequences, score guides, and rank candidates.

| API | Description |
|---|---|
| `design_guides` | Finds candidate guides for a target sequence and CRISPR system. |
| `find_pam_sites` | Scans a sequence for PAM matches. |
| `pam_matches` | Tests whether a sequence fragment satisfies a system PAM pattern. |
| `score_on_target` | Computes on-target activity score features. |
| `guide_gc_content` | Computes guide GC fraction. |
| `filter_guides` | Filters guide candidates by GC, scores, PAMs, or off-target summaries. |
| `rank_guides` | Orders guides by activity/specificity criteria. |

## 3. Off-Targets and Editing

Specificity and editing functions support practical guide selection.

| API | Description |
|---|---|
| `enumerate_off_targets` | Finds genomic off-target candidates up to mismatch limits. |
| `cfd_score` | Computes a CFD-like specificity score. |
| `mit_score` | Computes an MIT-like off-target score. |
| `guide_specificity_score` | Combines off-target candidates into guide-level specificity. |
| `genome_wide_offtarget_summary` | Summarizes off-target burden across genome sequences. |
| `design_base_editor_guides` | Designs guides whose editable bases fall inside an editor window. |
| `analyze_editing_window` | Reports target bases and bystander edits in an editing window. |
| `design_pegrna` | Designs prime-editing guide candidates. |
| `prime_editing_guide_score` | Scores pegRNA design features. |
| `design_hdr_template` | Builds HDR donor-template sequence suggestions. |
| `predict_indels` | Predicts likely indel outcomes around a cut site. |
| `indel_distribution` | Summarizes predicted indel lengths/classes. |

## 4. Screens and Libraries

Screen helpers analyze guide-level effects and library coverage.

| API | Description |
|---|---|
| `crispr_screen_analysis` | Aggregates guide counts/effects into gene-level screen results. |
| `mageck_like_test` | Runs a MAGeCK-like robust rank/count test. |
| `design_library` | Builds pooled guide libraries across targets. |
| `library_coverage_stats` | Computes target and guide coverage statistics for a library. |

---

## Complete Usage Example

```julia
using BioToolkit

system = SpCas9()
guides = design_guides(DNASeq("ACGT..."), system)
scored = rank_guides(guides)
offtargets = enumerate_off_targets(first(scored), genome_sequences)
```

