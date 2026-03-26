# metabolomics.jl

## Purpose
This file implements metabolomics-oriented analysis helpers on top of the proteomics layer. It focuses on source tracking, differential abundance, streaming analysis, and simple feature summaries for LC-MS-style data.

## Main type
- `MetabolomicsSourceTrackingResult` stores the raw chain plus mean, median, and credible-interval summaries for source proportions.

## Public functions
- `metabolomics_source_tracking_model(observed, source_profiles)` returns the Turing model used for source tracking.
- `metabolomics_source_tracking(observed, source_profiles; draws, rng)` runs sampling and returns posterior summaries.
- `metabolomics_source_tracking_posterior_summary(chain)` converts a chain into summary statistics.
- `metabolomics_differential_abundance(matrix, groups; covariates)` delegates to the proteomics differential-abundance path.
- `metabolomics_streaming_analysis(source, sink; threshold)` processes streamed experiments until a peak threshold is reached.
- `annotate_metabolite_features(features; labels)` returns a `DataFrame` with mean intensity and variance per feature.
- `quantify_metabolite_variation(features; robust=true)` returns a `DataFrame` with center and dispersion summaries.

## How it is used
This module is most useful when metabolomics data need source attribution or simple comparison across groups. A caller can feed observed counts and source profiles into `metabolomics_source_tracking`, use `metabolomics_differential_abundance` for group-level testing, and then summarize features with `annotate_metabolite_features` or `quantify_metabolite_variation`.

`metabolomics_streaming_analysis` is aimed at pipeline-style processing, where experiments are produced incrementally and downstream work should stop once a condition is met.

## Implementation notes
- The source-tracking routines require Turing; the module keeps the model and sampler references in lazy-loaded caches.
- The differential-abundance path reuses the proteomics implementation rather than duplicating model code.
- `quantify_metabolite_variation` can use either robust median/MAD-style summaries or standard mean/SD summaries.

## Why it matters
Metabolomics workflows often need Bayesian attribution and summary statistics that sit slightly outside generic proteomics. This file provides those workflows in a compact form while still reusing the shared mass-spectrometry infrastructure.
