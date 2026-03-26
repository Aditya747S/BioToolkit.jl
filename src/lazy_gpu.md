# lazy_gpu.jl

## Purpose
This file contains on-demand CUDA hooks for a few compute-heavy operations. It keeps GPU support optional by loading CUDA only when a GPU-backed path is actually needed.

## Internal state
- `_CUDA_SEQUENCE_LOADED`, `_CUDA_QUERY_LOADED`, and `_CUDA_PHYLO_LOADED` track whether the corresponding GPU specializations have been initialized.
- `_CUDA_SEQUENCE_HAMMING_IMPL`, `_CUDA_SEQUENCE_KMER_IMPL`, `_CUDA_SEQUENCE_DOTMATRIX_IMPL`, `_CUDA_QUERY_BIN_IMPL`, `_CUDA_QUERY_WINDOW_IMPL`, and `_CUDA_PHYLO_DISTANCE_MATRIX_IMPL` cache the lazily created implementations.

## Internal helpers
- `_is_cuda_backed_array` checks whether an array comes from CUDA.
- `_ensure_cuda_sequence!`, `_ensure_cuda_query!`, and `_ensure_cuda_phylo!` load CUDA and register specialized kernels on first use.

## What the GPU paths cover
- Sequence Hamming distance.
- GPU k-mer counting for short k values.
- Sequence dot-matrix construction.
- Query-position binning.
- Window-based coverage aggregation.
- Phylogenetic distance-matrix calculation for simple Hamming/p-distance use cases.

## How it is used
This file is not meant to be called directly by most users. Instead, other modules can check for CUDA-backed arrays and then call the lazy initializers to switch to the GPU implementation only when CUDA is present.

## Implementation notes
- CUDA is imported only inside the lazy initializer functions.
- The GPU kernels use atomic updates where multiple threads may write into the same histogram or coverage buffer.
- The phylogenetic distance matrix path falls back to a CPU-friendly zero matrix when the GPU branch is not applicable.

## Why it matters
The module lets BioToolkit keep a lightweight default install while still offering GPU acceleration for expensive routines. That balance is important in a package that needs to run both on small laptop setups and on larger accelerated systems.
