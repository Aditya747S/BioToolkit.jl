# `lazy_gpu.jl` - Lazy GPU and Threading Backend Helpers

## Overview

`lazy_gpu.jl` provides optional CUDA-backed implementations and lightweight backend helpers without making CUDA a hard dependency for BioToolkit's CPU-only workflows.

### Purpose

Many BioToolkit algorithms run well on CPU but can benefit from GPU acceleration for large sequence, query, and phylogenetic workloads. This file keeps GPU support lazy:

- CUDA-specific code is defined only when the CUDA extension has loaded.
- CPU workflows can load and run without CUDA installed.
- Backend resolution functions move data to GPU only when requested or available.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Lazy CUDA method creation** | CUDA kernels are created inside `_ensure_cuda_*` functions so CPU users do not pay CUDA load cost. |
| **Cached implementation refs** | `Ref{Any}` slots store CUDA implementation functions after first load. |
| **Backend detection by type/module name** | `_is_cuda_backed_array` avoids referencing CUDA types before CUDA is loaded. |
| **Small public helpers** | `cuda_available`, `resolve_backend`, `maybe_to_device`, and threading helpers provide a minimal backend abstraction. |
| **CPU fallback remains default** | `:auto` chooses CUDA only when it is actually available. |

---

## 1. CUDA Implementation Caches

```julia
const _CUDA_SEQUENCE_LOADED = Ref(false)
const _CUDA_SEQUENCE_HAMMING_IMPL = Ref{Any}(nothing)
const _CUDA_SEQUENCE_KMER_IMPL = Ref{Any}(nothing)
const _CUDA_SEQUENCE_DOTMATRIX_IMPL = Ref{Any}(nothing)

const _CUDA_QUERY_LOADED = Ref(false)
const _CUDA_QUERY_BIN_IMPL = Ref{Any}(nothing)
const _CUDA_QUERY_WINDOW_IMPL = Ref{Any}(nothing)

const _CUDA_PHYLO_LOADED = Ref(false)
const _CUDA_PHYLO_DISTANCE_MATRIX_IMPL = Ref{Any}(nothing)
```

**Kind:** Internal state

**Description:** These refs track whether CUDA implementations have been defined and cache function objects for later `invokelatest` calls.

---

## 2. CUDA Detection

### `_is_cuda_backed_array`

```julia
_is_cuda_backed_array(x) -> Bool
```

**Kind:** Internal predicate

**Description:** Returns true when `x` appears to be a CUDA `CuArray`, based on the type name and parent module name.

**Reason for implementation:** It avoids a direct hard reference to `CUDA.CuArray` when CUDA is not loaded.

---

## 3. Sequence CUDA Loader

### `_ensure_cuda_sequence!`

```julia
_ensure_cuda_sequence!() -> nothing
```

**Kind:** Internal lazy loader

**Description:** Defines and caches CUDA sequence helpers:

| Cached implementation | Purpose |
|---|---|
| `_CUDA_SEQUENCE_HAMMING_IMPL` | Hamming distance over `CuArray{UInt8,1}`. |
| `_CUDA_SEQUENCE_KMER_IMPL` | CUDA k-mer counting for DNA bytes. |
| `_CUDA_SEQUENCE_DOTMATRIX_IMPL` | CUDA dot-matrix match computation. |

**Important constraints:**

- Hamming distance requires equal sequence lengths.
- CUDA k-mer frequency requires `k > 0` and `k <= 10`.
- Ambiguous DNA bytes are skipped by the k-mer kernel.

---

## 4. Query CUDA Loader

### `_ensure_cuda_query!`

```julia
_ensure_cuda_query!() -> nothing
```

**Kind:** Internal lazy loader

**Description:** Defines and caches CUDA query helpers:

| Cached implementation | Purpose |
|---|---|
| `_CUDA_QUERY_BIN_IMPL` | Position binning into integer bins. |
| `_CUDA_QUERY_WINDOW_IMPL` | Per-window interval coverage. |

**Errors:** Window coverage requires `length(starts) == length(stops)`.

---

## 5. Phylogeny CUDA Loader

### `_ensure_cuda_phylo!`

```julia
_ensure_cuda_phylo!() -> nothing
```

**Kind:** Internal lazy loader

**Description:** Defines and caches a CUDA distance-matrix helper for phylogenetic workflows. Current GPU acceleration is focused on Hamming-like and p-distance calculations.

**Multi-device behavior:** When multiple CUDA devices and multiple Julia threads are available, sequence chunks are distributed across devices.

---

## 6. Backend Availability

### `cuda_available`

```julia
cuda_available() -> Bool
```

**Description:** Returns true if the CUDA extension is loaded and the runtime is functional.

**Use case:**

```julia
if cuda_available()
    # enable GPU path
end
```

### `flux_available`

```julia
flux_available() -> Bool
```

**Description:** Returns true when `Flux` is defined in the BioToolkit module, typically through an extension.

---

## 7. Backend Resolution

### `resolve_backend`

```julia
resolve_backend(; backend=:auto) -> Symbol
```

**Description:** Resolves a requested backend.

| Input | Result |
|---|---|
| `:auto` | `:cuda` if CUDA is available, otherwise `:cpu`. |
| `:cuda` | `:cuda`, but errors if CUDA is unavailable. |
| `:cpu` | `:cpu`. |

**Errors:** Unsupported backend symbols raise `ArgumentError`.

### `maybe_to_device`

```julia
maybe_to_device(x; backend=:auto)
```

**Description:** Moves `x` to a CUDA array when the resolved backend is `:cuda`. If `x` is already CUDA-backed, it is returned unchanged. CPU backend returns `x`.

### `maybe_to_host`

```julia
maybe_to_host(x)
```

**Description:** Converts CUDA-backed arrays to host `Array`; returns CPU values unchanged.

---

## 8. Threading Helpers

### `threaded_foreach`

```julia
threaded_foreach(count::Integer, f::Function; threaded=true)
```

**Description:** Applies `f(i)` for `i in 1:count`, using `Threads.@threads` when requested and more than one thread is available.

### `threaded_map_collect`

```julia
threaded_map_collect(f::Function, xs; threaded=true)
```

**Description:** Collects `f(x)` for each item in `xs`, optionally using threaded execution.

---

## Quick Reference

| API | Purpose |
|---|---|
| `_ensure_cuda_sequence!` | Define sequence CUDA kernels lazily. |
| `_ensure_cuda_query!` | Define query CUDA kernels lazily. |
| `_ensure_cuda_phylo!` | Define phylogeny CUDA helper lazily. |
| `cuda_available` | Check CUDA runtime availability. |
| `flux_available` | Check Flux extension availability. |
| `resolve_backend` | Resolve `:auto`, `:cpu`, or `:cuda`. |
| `maybe_to_device` | Move data to CUDA when requested. |
| `maybe_to_host` | Convert CUDA-backed data to CPU. |
| `threaded_foreach` | Optional threaded loop. |
| `threaded_map_collect` | Optional threaded map/collect. |

---

## Complete Usage Example

```julia
backend = resolve_backend(backend=:auto)
x_dev = maybe_to_device([1, 2, 3]; backend=backend)
x_host = maybe_to_host(x_dev)

values = threaded_map_collect(1:10; threaded=true) do x
    x^2
end
```
