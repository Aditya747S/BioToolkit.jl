# `coevolution.jl` - Coevolution and Contact Prediction

## Overview

`coevolution.jl` analyzes multiple sequence alignments for conservation, covariation, direct coupling, and residue contact prediction.

### Purpose

This page documents the user-facing types and workflows implemented in `coevolution.jl`. The emphasis is on public APIs exported through BioToolkit and the biological analysis step each one supports.

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

## 1. Core Types

These structures hold contact maps and fitted coupling models.

| API | Description |
|---|---|
| `ContactMap` | Predicted or observed residue contact pairs with scores and optional metadata. |
| `PseudoLikelihoodModel` | Model parameters from pseudo-likelihood direct coupling analysis. |

## 2. DCA Workflow

A complete contact-prediction workflow starts with alignment filtering and ends with ranked residue pairs.

| API | Description |
|---|---|
| `filter_alignment_for_dca` | Removes poorly covered sequences/columns before coupling analysis. |
| `sequence_reweighting` | Computes sequence weights to reduce redundancy in deep alignments. |
| `fit_pseudolikelihood_model` | Fits a pseudo-likelihood model for direct coupling analysis. |
| `compute_contact_scores` | Converts fitted couplings into residue-pair contact scores. |
| `predict_contact_map` | Builds a `ContactMap` from alignment/model scores. |
| `top_contact_pairs` | Returns the highest-scoring nonlocal residue pairs. |
| `fold_from_contacts` | Creates a coarse fold/contact graph from predicted contacts. |

## 3. Covariation and Diagnostics

Diagnostic functions expose alternative scores and alignment quality signals.

| API | Description |
|---|---|
| `mutual_information_contacts` | Scores residue pairs by mutual information. |
| `direct_information_contacts` | Computes direct-information-style coupling scores. |
| `column_conservation_scores` | Computes per-column conservation. |
| `sequence_logo_entropy` | Computes entropy values suitable for sequence logos. |
| `evolutionary_coupling_network` | Builds a network representation of strong couplings. |
| `contact_enrichment_statistics` | Summarizes contact enrichment against a supplied reference. |
| `shrinkage_precision_contacts` | Estimates contacts from shrinkage precision/covariance. |
| `phylogenetic_correction` | Applies correction for phylogenetic relatedness. |
| `positional_covariation_matrix` | Builds a full position-by-position covariation matrix. |
| `gap_analysis` | Reports gap burden by sequence and column. |
| `contact_precision_recall` | Compares predicted contacts against known contacts. |
| `alignment_quality_report` | Combines depth, gaps, diversity, and conservation diagnostics. |

---

## Complete Usage Example

```julia
using BioToolkit

filtered = filter_alignment_for_dca(msa)
weights = sequence_reweighting(filtered)
model = fit_pseudolikelihood_model(filtered; weights=weights)
contacts = predict_contact_map(filtered, model)
top = top_contact_pairs(contacts, 20)
```

