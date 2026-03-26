# Examples

This folder contains small runnable examples for the repository-level `BioToolkit` package.

## What to run first

1. `dna_workflow.jl`
2. `vcf_workflow.jl`
3. `bed_workflow.jl`
4. `gpu_ready_structs.jl`
5. `restriction_workflow.jl`
6. `entrez_workflow.jl`
7. `kegg_pathway_workflow.jl`
8. `Comparison/cuda_histogram_compare.jl`

Each `.jl` file has a matching `.md` file with a detailed explanation written for readers who may be new to Julia.

The benchmark-focused files live under `Examples/Comparison/`; see that folder's README for a short guide.

## How to run an example

From the repository root, run:

```bash
julia --project=. Examples/vcf_workflow.jl
```

Replace the file name with the example you want to try.

## What these examples cover

- parsing and ingesting VCF text
- DNA counting, transcription, translation, ORFs, and k-mers
- loading Arrow output back into Julia
- filtering and histogramming genomic coordinates
- ingesting BED-like interval data
- preparing data that is compatible with GPU transfer when CUDA.jl is available
- restriction-enzyme scanning and simple in-memory digest workflows
- live Entrez / NCBI search and fetch examples
- KEGG pathway graph extraction and Mermaid export
- comparing Julia CPU, Julia CUDA, and Python CPU histogram performance on synthetic genomic positions
