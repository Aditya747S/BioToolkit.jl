# `Examples/restriction_workflow.jl`

This example shows the restriction-enzyme workflow in BioToolkit:

1. inspect the enzyme catalog
2. look up a named enzyme
3. scan a DNA sequence for cut sites
4. digest a sequence with multiple enzymes
5. inspect the enzyme-to-hit mapping helper

## Why this example matters

Restriction analysis is one of the most common wet-lab bioinformatics tasks. The package now includes a built-in catalog of common enzymes, search helpers for finding cut sites, and a digestion helper for simple in-memory workflows.

## What the code demonstrates

### `restriction_enzyme_names`

Returns the catalog names in sorted order.

### `restriction_enzyme`

Looks up a named enzyme and returns a compact `RestrictionEnzyme` value.

### `find_restriction_sites`

Scans a DNA sequence for matches to a specific enzyme recognition site.

### `digest_sequence`

Produces fragments after simulating one or more cuts.

### `restriction_digest_map`

Returns a per-enzyme hit map, which is useful when you want to inspect several enzymes at once without recomputing your own loops.

## How to run it

From the repository root:

```bash
julia --project=. Examples/restriction_workflow.jl
```

## Notes

The built-in catalog focuses on common lab enzymes first. The module is intentionally lightweight and keeps the matching logic in-memory, so it stays fast for short primer-scale and plasmid-scale workflows.
