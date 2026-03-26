# `Examples/entrez_workflow.jl`

This example shows the live NCBI / Entrez workflow in BioToolkit:

1. search PubMed
2. summarize returned hits
3. search taxonomy
4. search nucleotide records
5. fetch FASTA and GenBank text for a returned accession

## Why this example matters

Biopython users often rely on `Bio.Entrez` as the first step in a workflow: search the public databases, inspect a few records, and then fetch the sequence or annotation payload you actually need. BioToolkit now includes the same style of workflow with compact helper functions.

## What the code demonstrates

### `entrez_pubmed_search`

Searches PubMed using the public E-utilities endpoint.

### `entrez_summary`

Fetches compact JSON summaries for a set of IDs.

### `entrez_taxonomy_search`

Convenience wrapper for the taxonomy database.

### `entrez_nucleotide_search`

Searches NCBI nucleotide / nuccore records.

### `entrez_fetch_fasta`

Fetches a FASTA payload for a returned accession.

### `entrez_fetch_genbank`

Fetches the GenBank flat-file payload for a returned accession.

## How to run it

From the repository root:

```bash
julia --project=. Examples/entrez_workflow.jl
```

## Notes

- This is a live network example, so it depends on NCBI being reachable.
- The example is wrapped in `try` / `catch` so it prints a useful message rather than crashing if the network is unavailable or NCBI rate limits the request.
- If you use this in real workflows, set a descriptive email and tool name through the keyword arguments on the Entrez helpers.
