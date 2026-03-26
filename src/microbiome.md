# microbiome.jl

## Purpose
This file implements microbial community analysis on top of count matrices, phylogenetic trees, and sample metadata. It covers compositional transforms, beta-diversity, ordination, differential abundance, source tracking, and network construction.

## Main types
- `CommunityProfile` is the central container combining counts, taxonomy, a phylogenetic tree, and metadata.
- `PCoAResult` stores ordination coordinates, eigenvalues, and explained variance.
- `NMDSResult` stores ordination coordinates, stress, iterations, and convergence state.
- `ANCOMResult` stores taxon IDs, W statistics, p-values, q-values, and significance flags.
- `SongbirdResult` stores regression coefficients and fit metadata for rank-based differential abundance.
- `SourceTrackingResult` stores posterior summaries from a source tracking model.
- `MicrobiomeNetwork` stores a graph, edge weights, taxon labels, and coordinates.

## Key functions
- Compositional transforms: `clr_transform` and `ilr_transform`.
- Distance and diversity metrics: `bray_curtis`, `pairwise_bray_curtis`, `unifrac`, `pairwise_unifrac`, `weighted_unifrac`, `shannon_entropy`, `simpson_index`, and `faith_pd`.
- Ordination: `pcoa`, `pcoa_plot`, and `nmds`.
- Differential abundance: `ancom` and `songbird`.
- Source tracking: `source_tracking_model`, `source_tracking`, and `source_tracking_posterior_summary`.
- Networks: `cooccurrence_network` and `network_plot`.

## How it is used
The normal workflow is to build a `CommunityProfile`, transform counts with `clr_transform` or `ilr_transform`, compute distances with `bray_curtis` or `unifrac`, and then reduce the result with `pcoa` or `nmds`.

For hypothesis testing, `ancom` and `songbird` provide two different differential-abundance style approaches. For source attribution, the source tracking helpers run a Bayesian model and summarize the posterior chain. Network functions turn taxa associations into graph objects that can be plotted or inspected downstream.

## Implementation notes
- The module validates that taxonomy and metadata line up with the count matrix and phylogeny before constructing a `CommunityProfile`.
- `pairwise_unifrac` supports both weighted and unweighted workflows.
- `pcoa_plot` expects at least two coordinates and can annotate points with labels.
- `integrate_batches`, `detect_doublets`, and clustering support are implemented in the downstream single-cell module, not here.

## Why it matters
Microbiome analysis is not just counting taxa; it is compositional, phylogenetic, and often posterior-driven. This file pulls those pieces into one coherent API so community analysis is reusable instead of being reimplemented in ad hoc scripts.
