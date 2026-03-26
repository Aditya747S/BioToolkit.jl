# BioToolkit.jl

> High-performance bioinformatics toolkit for Julia with a unified API for sequence analysis, structural biology, phylogenetics, genomics, and omics workflows.

[Benchmark report](benchmark_report.md) | [API reference](api.md)

[![License: Custom](https://img.shields.io/badge/License-Custom-yellow.svg)](https://github.com/Aditya747S/BioToolkit.jl/blob/main/LICENSE)
[![Julia Version](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org/)

BioToolkit features original, optimized implementations of core algorithms. It is designed as a cohesive toolkit rather than a bundle of wrappers, so sequence I/O, phylogenetics, structural biology, and omics workflows all live in one namespace.

## Highlights

| Area | What you get |
|:--|:--|
| Sequence and protein analysis | GC content, translation, codon usage, melting temperature, quality scores, and fast k-mer workflows. |
| Structural biology | PDB and mmCIF parsing, RMSD, Kabsch superposition, SASA, contact maps, and residue-level geometry. |
| Phylogenetics and PopGen | Neighbor-joining, consensus trees, phylogenetic distances, $F$-statistics, GenePop, and GWAS scanning. |
| Omics and systems biology | Differential expression, single-cell analysis, epigenetics, metabolomics, microbiome, and clinical workflows. |

## Performance snapshot

| Workflow | BioToolkit | Baseline | Speedup |
| :--- | :--- | :--- | :--- |
| Hamming distance (19.5 kbp) | 0.79 ms | Biopython: 259.6 ms | 328x |
| K-mer frequency (CUDA) | 0.14 ms | Biopython: 302.6 ms | 2161x |
| BED parsing | 7.7 ms | BED.jl: 38.2 ms | 5.0x |
| GenBank parsing (500 records) | 3.8 ms | Biopython: 21.4 ms | 5.7x |
| Phylogenetics (NJ tree, 200 taxa) | 11.0 ms | Biopython: 1848.5 ms | 168x |
| Structural superposition | 0.06 ms | Bio.PDB: 0.75 ms | 12.6x |

See the full [benchmark report](benchmark_report.md) for details.

## Quick start

```julia
using BioToolkit

seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

gc_content(seq)
translate_dna(seq)
reverse_complement(seq)
kmer_frequency(seq, 3)
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Aditya747S/BioToolkit.jl.git")
```

## Optional backends

* `BioToolkitPlotsExt` for Plots.jl
* `BioToolkitMakieExt` for Makie
* `BioToolkitTuringExt` for probabilistic modeling hooks

## Project status

BioToolkit is in active development. The package is built as a single coherent toolkit, so the docs focus on integrated workflows instead of isolated wrappers.

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request with tests and updated documentation.

## License

Please see the [LICENSE](https://github.com/Aditya747S/BioToolkit.jl/blob/main/LICENSE) file for details.