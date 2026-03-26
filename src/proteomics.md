# proteomics.jl

## Purpose
This file implements a compact mass-spectrometry proteomics workflow. It starts with mzML parsing and ends with peak detection, alignment, imputation, statistical testing, sparse classification, streaming, and de novo sequence reconstruction.

## Main types
- `Spectrum` stores one scan with retention time, m/z values, and intensities.
- `MassSpecExperiment` bundles spectra, a feature matrix, and a sample design table.
- `MassSpecPeak` stores a detected chromatographic or spectral peak.
- `PeakDetectionResult` stores the peaks together with multiscale coefficients and the scales used.
- `AlignmentResult` stores the warping path, warped query signal, total cost, and warp vector.
- `DifferentialAbundanceResult` stores coefficients, test statistics, q-values, and feature/sample labels.
- `SparsePLSDAResult` stores weights, scores, loadings, VIP scores, selected feature indices, and class labels.
- `DeNovoResult` stores the reconstructed peptide sequence, graph path, score, and spectrum graph.

## Supporting helpers
- `_AA_MASSES` maps amino-acid letters to monoisotopic masses for spectrum graph construction.
- `_mad`, `_parse_float_vector`, `_extract_tag_value`, `_extract_series`, and `_parse_spectrum` are internal parsing and robust-statistics helpers.
- `_ricker` and `_same_convolution` implement the multiscale wavelet-style peak detector.
- `_dtw_alignment`, `interp1d`, and `_design_matrix` provide alignment and model-building support.
- `qrilc_impute` uses a log-normal style truncated distribution to fill missing or nonpositive values.

## Public functions
- `read_mzml(path)` parses a simplified mzML file into a `MassSpecExperiment`.
- `detect_peaks(signal; scales, threshold, mz, rt, use_gpu)` detects peaks from a signal or spectrum.
- `align_samples(reference, query; method=:dtw, band=0)` aligns two sample traces using DTW or a simple obiwarp-like transform.
- `qrilc_impute(matrix; rng)` imputes missing proteomics values.
- `differential_abundance(matrix, groups; covariates, qvalue_threshold)` fits a per-feature linear model and adjusts p-values.
- `mixed_model_abundance(matrix, groups; batch=nothing)` currently delegates to differential abundance.
- `sparse_pls_da(X, y; n_components, sparsity)` performs sparse PLS-DA-style feature selection.
- `stream_mass_spec(producer, consumer; chunk_size)` wires a producer and consumer together through a channel.
- `build_spectrum_graph(masses; tolerance)` converts mass gaps into a directed graph.
- `de_novo_sequence(masses; tolerance)` reconstructs the best peptide path from a spectrum graph.

## How it is used
The usual workflow is: parse data with `read_mzml`, detect features with `detect_peaks`, align traces with `align_samples`, impute missing intensities with `qrilc_impute`, then run `differential_abundance` or `sparse_pls_da` on the resulting feature matrix.

`build_spectrum_graph` and `de_novo_sequence` support a separate de novo interpretation path for tandem-mass-spectrum style data. `stream_mass_spec` is useful when data should be processed incrementally instead of loaded eagerly.

## Implementation notes
- The module depends on `DifferentialExpression.benjamini_hochberg` for q-value adjustment.
- `read_mzml` is intentionally lightweight and uses line-oriented parsing rather than a full XML stack.
- `detect_peaks` uses a robust median and MAD cutoff, which is stable on noisy spectra.
- `mixed_model_abundance` is a compatibility wrapper rather than a full mixed-model fit.

## Why it matters
This module gives BioToolkit a complete proteomics analysis path in one place. It connects raw spectra to downstream statistical and classification workflows instead of leaving those steps as separate scripts.
