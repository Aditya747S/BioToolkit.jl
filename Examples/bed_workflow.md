# `Examples/bed_workflow.jl`

This example shows how to work with BED-like interval data.

## What this example covers

1. write a small BED file
2. ingest it into Arrow
3. load the Arrow file back into Julia
4. inspect the interval columns
5. check how chunked ingestion affects the Arrow file

## Who should read this

This example is useful if you are new to BED files, new to Julia, or simply want the shortest path from a plain text interval file to a fast columnar representation.

## What the code demonstrates

### `parse_bed_record`

The package parses each BED line into a `BedRecord`.

A BED record in this package is intentionally minimal:

- `chrom`
- `start`
- `stop`

That keeps the code easy to understand and easy to extend.

### `ingest_bed`

`ingest_bed` reads the BED file in chunks and writes Arrow record batches.

The point is the same as with VCF:

- text parsing happens once
- Arrow is used for repeated analysis

### `load_arrow_table`

After ingestion, the Arrow file can be loaded and inspected as a columnar table.

The example calculates the total interval length directly from the loaded columns:

```julia
table.stop .- table.start
```

That works because the Arrow table gives you column vectors, and Julia can do simple arithmetic on those vectors directly.

To count rows, use `Tables.rowcount(table)` or `length(table.start)`. This is better than `length(table)` because table objects often report columns when you ask for their length.

### `Arrow.Stream`

The example checks the number of partitions written to Arrow.

This is a nice sanity check for chunked ingestion because it shows that the file was written in multiple batches instead of as one large batch.

### `window_coverage`

The example also computes fixed-width window coverage for `chr1`.

This is the smallest useful summary if you want to reproduce the kind of windowed coverage plots that show up in many genomics papers:

- each window has a fixed width
- each interval contributes only the bases that overlap that window
- the output is easy to inspect or plot later

For the sample file in this example, the windows on `chr1` are all 10 bases wide, so the result is a simple check that the interval arithmetic is working correctly.

## How to run it

From the repository root:

```bash
julia --project=. Examples/bed_workflow.jl
```

## How to modify it

To adapt the example to your own BED file:

- replace the sample writer with your real input path
- keep `ingest_bed` the same
- use `load_arrow_table` to inspect the result
- use the loaded columns for your own interval calculations

## Notes on GPU usage

This example is primarily CPU-side because BED data contains strings in the `chrom` column.

That means:

- `BedRecord` is useful for ingestion and file conversion
- `BedRecord` is not the GPU-ready event type
- the numeric interval columns `start` and `stop` are easy to work with on the CPU

If you later want a GPU workflow, the usual pattern is to encode chromosomes into numeric IDs before sending the data to device memory.
