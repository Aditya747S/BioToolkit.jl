# `Examples/gpu_ready_structs.jl`

This example explains which package types are ready for GPU-style usage and why.

## The main idea

The package separates two kinds of data:

- parsing records that are easy to read and easy to ingest
- compact event records that are suitable for high-performance computation

That distinction matters if you want to use CUDA.jl later.

## What the code shows

### `VariantTextRecord`

`VariantTextRecord` is the CPU-side parsing record for VCF text.

It contains strings, so it is convenient for text processing but not ideal for GPU memory.

### `VariantEvent`

`VariantEvent` is the compact event type.

It uses small numeric fields instead of strings, so it is an isbits type. That makes it a much better candidate for GPU transfer.

### `compact_variant_event`

This function converts a parsed text record into the compact numeric form.

That is the key step if you want to move from file parsing to fast computation.

### `BedRecord`

`BedRecord` is also CPU-side because it stores a chromosome string.

That is perfectly fine for ingestion and file conversion, but it is not the final GPU-friendly type.

## How to run it

From the repository root:

```bash
julia --project=. Examples/gpu_ready_structs.jl
```

## How to use CUDA with these types

The repository itself does not depend on CUDA.jl. That is intentional.

If you install CUDA.jl in your own environment, the usual workflow is:

1. parse and compact on the CPU
2. create a `Vector{VariantEvent}`
3. copy that vector to the GPU
4. run your GPU code on the device array

Example pattern:

```julia
using CUDA
using BioToolkit

events = BioToolkit.compact_variant_event.([
    BioToolkit.parse_vcf_record("chr1\t10\trs1\tA\tG\t99.5"),
])

gpu_events = CuArray(events)
```

Why this works:

- `VariantEvent` is isbits-compatible
- GPU memory transfer works best with fixed-size, concrete data
- strings and other heap-managed objects are avoided in the device-side representation

## Important caution

Only the compact event type is the GPU-friendly one.

These are not GPU-ready in their current form:

- `VariantTextRecord`
- `BedRecord`

They are still useful, but only for CPU-side parsing and ingestion.

## What a student should take away

If you are new to Julia or GPUs, the simplest mental model is this:

- use the text record types to read files
- convert to `VariantEvent` when you need speed
- move `Vector{VariantEvent}` to CUDA.jl if you want GPU computation

That keeps the pipeline understandable and minimizes the amount of code you need to change.
