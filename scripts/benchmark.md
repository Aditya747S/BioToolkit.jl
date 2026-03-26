# `scripts/benchmark.jl`

This file is a tiny benchmark harness for the current fast paths in the repository.

## What it does

`benchmark.jl`:

- generates small synthetic VCF and BED inputs
- times chunked ingestion into Arrow
- times region filtering
- times histogram binning in both serial and threaded modes
- prints the results in milliseconds

## Why it matters

The repository is focused on speed, so having a small benchmark script is useful. It lets you check the current performance of the code without introducing a heavy benchmarking framework.

## Design goals

The script is intentionally minimal:

- standard library only, plus the package itself
- no extra abstractions
- easy to read and modify

That makes it suitable for students who want to change the sample size, alter the interval ranges, or compare their own changes against the current baseline.

## What it measures

The script times the current core paths:

- `ingest_vcf`
- `ingest_bed`
- `filter_region`
- `bin_positions` with `threaded=true`
- `bin_positions` with `threaded=false`

For DNA sequence work, see `scripts/sequence_benchmark.jl`, which compares the Julia toolkit against Biopython.

For native Julia package comparisons, see `Examples/Comparison/julia_package_compare.jl`, which covers BioSequences, BioAlignments, BED, and GFF3 workloads.
