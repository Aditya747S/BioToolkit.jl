# `deeplearning.jl` - Deep Learning Feature Helpers

## Overview

`deeplearning.jl` provides model-inspired embedding, denoising, label-transfer, perturbation, protein-sequence, and graph-network utilities using lightweight Julia implementations or payload-compatible approximations.

### Purpose

This page documents the user-facing types and workflows implemented in `deeplearning.jl`. The emphasis is on public APIs exported through BioToolkit and the biological analysis step each one supports.

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

## 1. Single-Cell Embeddings

Embedding helpers produce latent spaces for downstream clustering, integration, or visualization.

| API | Description |
|---|---|
| `scvi_like_embedding` | Builds a variational-style low-dimensional embedding from count data. |
| `geneformer_like_embedding` | Computes gene-token-inspired cell embeddings. |
| `scgpt_like_embedding` | Builds transformer-inspired gene-expression embeddings. |
| `scbert_like_embedding` | Builds BERT-style expression embeddings. |
| `batch_corrected_latent` | Produces latent features adjusted for batch labels. |
| `contrastive_cell_embedding` | Learns or approximates contrastive cell embeddings. |

## 2. Prediction and Transfer

These functions map labels, classify cells, denoise data, or predict perturbation responses.

| API | Description |
|---|---|
| `cellassign_like_mapping` | Maps cells to marker-defined cell types. |
| `graphsca_label_transfer` | Transfers labels over a graph/embedding neighborhood. |
| `zero_shot_cell_annotation` | Assigns cell types from marker descriptions without fitted labels. |
| `flux_mlp_classifier` | Trains or applies a compact Flux MLP classifier. |
| `cell_type_denoising` | Denoises expression by cell type or neighborhood. |
| `scgen_like_perturbation` | Predicts perturbation shifts in latent/expression space. |

## 3. Representations and Networks

General representation helpers support multimodal, protein, and GRN workflows.

| API | Description |
|---|---|
| `flux_autoencoder_embedding` | Learns an autoencoder latent representation. |
| `sparse_autoencoder_features` | Extracts sparse autoencoder feature activations. |
| `trajectory_neural_ode` | Fits or summarizes neural-ODE-like trajectory dynamics. |
| `multimodal_wnn_embedding` | Combines modalities in a WNN-style latent space. |
| `protein_sequence_embedding` | Computes protein sequence embeddings/features. |
| `attention_grn` | Infers regulatory edges from attention-like feature weights. |
| `gene_regulatory_network_gnn` | Builds graph-neural-network-inspired GRN scores. |
| `self_supervised_pretraining` | Runs self-supervised pretraining utilities for expression matrices. |
| `deep_factorization_embedding` | Computes nonlinear factorization embeddings. |
| `cell_cycle_regression` | Regresses or scores cell-cycle effects. |

---

## Complete Usage Example

```julia
using BioToolkit

latent = scvi_like_embedding(counts; n_latent=16)
labels = zero_shot_cell_annotation(counts, marker_sets)
perturbed = scgen_like_perturbation(latent, condition_labels)
```

