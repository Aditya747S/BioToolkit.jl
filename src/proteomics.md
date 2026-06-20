# `proteomics.jl` - Proteomics and Mass Spectrometry

## Overview

`proteomics.jl` supports mass-spectrometry experiment containers, spectra and peaks, mzML-style reading, peak detection, sample alignment, missing-value imputation, differential abundance, mixed models, sparse PLS-DA, DIA quantification, PTM localization/enrichment, peptide-to-structure projection, protein inference, spectrum graphs, and de novo sequencing.

### Purpose

This page documents the public BioToolkit APIs implemented in `proteomics.jl`. The sections follow the main workflow stages and describe how each function or type is intended to be used.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Public API coverage** | Exported types and functions are grouped by use case. |
| **Concrete descriptions** | Entries describe real behavior rather than generic call placeholders. |
| **Analysis-ready outputs** | Results are suitable for downstream visualization, statistics, or reporting. |
| **Interoperable inputs** | APIs work with standard Julia matrices/tables and BioToolkit containers. |
| **Composable workflows** | Higher-level functions build on lower-level summaries where possible. |

---

## 1. Core Types

Proteomics data are represented by experiment, spectrum, peak, and result containers.

| API | Description |
|---|---|
| `MassSpecExperiment` | Experiment-level container for spectra, sample metadata, and feature tables. |
| `Spectrum` | One MS/MS or MS1 spectrum with m/z, intensity, precursor, and metadata. |
| `MassSpecPeak` | Detected chromatographic or spectral peak. |
| `PeakDetectionResult` | Peak list plus detection parameters and diagnostics. |
| `AlignmentResult` | Retention-time/sample alignment result. |
| `DifferentialAbundanceResult` | Feature-level differential abundance result. |
| `SparsePLSDAResult` | Sparse PLS-DA model/result. |
| `DeNovoResult` | De novo peptide sequencing result. |

## 2. I/O and Preprocessing

These functions read spectra, detect peaks, align samples, and handle missingness.

| API | Description |
|---|---|
| `read_mzml` | Reads mzML-like mass spectrometry data. |
| `stream_mass_spec` | Streams spectra/features from large files. |
| `detect_peaks` | Detects peaks from spectra or chromatograms. |
| `align_samples` | Aligns peaks/features across samples. |
| `qrilc_impute` | Imputes missing left-censored intensity values with QRILC-style behavior. |

## 3. Statistics and Advanced Workflows

Analysis helpers cover abundance modeling, classification, DIA, PTMs, and peptide interpretation.

| API | Description |
|---|---|
| `differential_abundance` | Tests protein/peptide/feature differential abundance. |
| `mixed_model_abundance` | Fits mixed models for abundance with random/fixed effects. |
| `sparse_pls_da` | Runs sparse PLS-DA for class separation. |
| `build_spectrum_graph` | Builds graph representation of spectrum peaks/transitions. |
| `de_novo_sequence` | Infers peptide sequences from MS/MS spectra. |
| `dia_like_quantification` | Quantifies peptides/proteins from DIA-like data. |
| `phosphosite_localization` | Scores phosphosite localization. |
| `glycoproteomics_motif_table` | Builds glycosylation motif summaries. |
| `project_peptides_to_structure` | Maps peptide evidence onto protein structures. |
| `protein_inference_top3` | Infers protein abundance by top-three peptides. |
| `ptm_site_enrichment` | Tests PTM site enrichment across groups or annotations. |

---

## Complete Usage Example

```julia
using BioToolkit

exp = read_mzml("run.mzML")
peaks = detect_peaks(exp)
aligned = align_samples([peaks])
da = differential_abundance(feature_matrix, group_labels)
```

