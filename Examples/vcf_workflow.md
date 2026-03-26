# `Examples/vcf_workflow.jl`

This example shows the most important end-to-end workflow in the package:

1. create a small VCF file
2. ingest it into Arrow
3. load the Arrow file back into Julia
4. filter a genomic region
5. compute a histogram
6. inspect the record-batch partitioning written by Arrow

## Who this is for

This file is written for readers who may be new to Julia and new to the package. It avoids assuming prior knowledge of the Julia package manager, Arrow, or genomic file formats.

## What the code demonstrates

### `using BioToolkit`

This loads the package from the repository-level Julia project.

### `ingest_vcf`

`ingest_vcf` reads a VCF-like text file and writes it to Arrow in chunks.

Why that matters:

- parsing text is slow
- Arrow is columnar and much faster to reuse
- chunked writes keep the code memory-friendly

The example uses a small sample file so you can run it immediately, but the same code path is the one you would use for larger real data.

### `load_arrow_table`

`load_arrow_table` opens the Arrow file so you can read it back without reparsing the original text.

This is the main speed win in the package:

- ingest once
- analyze many times

When you want the number of rows in a table, use `Tables.rowcount(table)` or count the length of a specific column such as `length(table.pos)`. In Julia, `length(table)` can mean the number of columns for some table types, which is not what beginners usually expect.

### `filter_region`

`filter_region` returns rows in a chromosome and coordinate range.

The example shows both modes:

- the normal mode, which scans the table
- the sorted mode, which uses binary search when the input is already ordered

Use `sorted=true` only when the data is actually sorted by chromosome and position.

### `coverage_histogram`

`coverage_histogram` groups positions into bins.

This is the simplest possible coverage-style analysis and is intentionally easy to read.

### `Arrow.Stream`

The example also checks how many Arrow partitions were written. That shows the chunked writer is producing multiple record batches, which is what you want when you are trying to keep ingestion memory-friendly.

## How to run it

From the repository root:

```bash
julia --project=. Examples/vcf_workflow.jl
```

## How to modify it

If you want to adapt this example to your own VCF file:

- replace the sample writer with your input path
- keep `ingest_vcf(input_path, output_path)` the same
- point `load_arrow_table` at the Arrow file you created
- use `filter_region` for chromosome window queries
- use `coverage_histogram` for simple bin counts

## GPU notes

The VCF example itself runs entirely on the CPU, but it produces the data shape that can later be moved to a GPU.

The important type is `VariantEvent`:

- it is compact
- it uses only isbits-compatible fields
- it is designed so a vector of events can be copied into GPU memory if CUDA.jl is available

In contrast:

- `VariantTextRecord` is a parsing type and contains `String` fields
- `BedRecord` is also CPU-side because it stores strings

Those two are useful on the CPU, but `VariantEvent` is the type you would hand off to GPU code.
