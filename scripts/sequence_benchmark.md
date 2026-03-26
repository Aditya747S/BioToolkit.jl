# `scripts/sequence_benchmark.jl`

This script compares the small Julia DNA toolkit against Biopython.

## What it measures

The benchmark times the same kinds of operations from the playlist-style DNA toolkit:

- GC content
- reverse complement
- translation

It runs them on a repeated synthetic DNA sequence so the work is large enough to measure.

## Why Biopython

Biopython is the most direct Python baseline for this kind of sequence work.

It has mature sequence utilities, and it is the clearest comparison point for the functions added in this repository.

## Why the benchmark is structured this way

The script does two things:

1. warms up the Julia functions first
2. repeats each operation many times before timing

That keeps the comparison focused on actual execution cost rather than first-call compilation noise.

## What to expect

The Julia implementation should be competitive or faster on these small deterministic sequence transforms, especially when the same string work is repeated many times.

## How to run it

From the repository root:

```bash
julia --project=. scripts/sequence_benchmark.jl
```

The script also runs a short Python benchmark through the configured Conda environment so the numbers appear side by side.