# `Examples/kegg_pathway_workflow.jl`

This example shows the KEGG pathway workflow in BioToolkit:

1. parse a KEGG-style pathway flat file
2. build a graph representation
3. render the graph as Mermaid
4. write the Mermaid output to a file
5. compare the generated output against a checked-in fixture

## Why this example matters

KEGG-style pathway files are a common source of biological context, and a compact graph export makes them easy to inspect or embed in docs. The example stays offline and uses a small fixture so the output remains deterministic.

## What the code demonstrates

### `read_kegg_pathway`

Parses the pathway record into a structured `KEGGPathwayRecord`.

### `read_pathway_graph`

Converts the record into a graph with pathway, gene, enzyme, and compound nodes.

### `kegg_pathway_mermaid`

Renders the graph to Mermaid text.

### `write_kegg_pathway_mermaid`

Writes Mermaid output to disk, which is useful if you want to commit a diagram or inspect it in a browser-based Mermaid renderer.

## How to run it

From the repository root:

```bash
julia --project=. Examples/kegg_pathway_workflow.jl
```

## Notes

- The example compares the output against `Examples/fixtures/kegg_pathway_expected.mmd`.
- This is intentionally small and deterministic so it can act as a regression check for the renderer.
