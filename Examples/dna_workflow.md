# `Examples/dna_workflow.jl`

This example mirrors the early playlist topics: counting nucleotides, transcription, reverse complements, GC content, translation, ORFs, and k-mers.

## What the code demonstrates

### `validate_dna`

This checks that a sequence uses only DNA letters supported by the toolkit.

That is the smallest useful guard you want before calling the rest of the sequence functions.

### `count_nucleotides`

This returns a named tuple with counts for `A`, `C`, `G`, `T`, and `N`.

It is the direct answer to the classic beginner bioinformatics task of counting bases in a DNA string.

### `transcribe_dna`

This converts DNA to RNA by replacing `T` with `U`.

That matches the standard transcription exercise from the playlist.

### `reverse_complement`

This computes the reverse complement.

The implementation stays simple and direct so the control flow is easy to follow.

### `gc_content`

This returns the GC fraction as a number between 0 and 1.

That is the same summary used in many beginner exercises and in real sequence QC.

### `translate_dna`

This walks through the sequence in codons and translates them into amino acids.

It is the core of the translation section in the playlist.

### `find_orfs`

This produces simple ORF candidates by translating each frame.

It is intentionally lightweight rather than a full gene finder.

### `kmer_frequency`

This counts fixed-length subsequences.

That is the same pattern behind basic genome statistics and pattern discovery.

## How to run it

From the repository root:

```bash
julia --project=. Examples/dna_workflow.jl
```

## Why this is fast enough to matter

The code is small, concrete, and allocation-light:

- no dataframes
- no Python object overhead
- no extra parsing framework
- simple loops over strings

That is the kind of structure that Julia can optimize well.

## How to extend it

If you want to grow this example, the next natural additions are:

- reading a real FASTA file with `read_fasta`
- computing k-mers for multiple values of `k`
- searching for ORFs in reverse-complement frames
- wrapping the outputs in plots or tabular summaries