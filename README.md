# BioToolkit.jl

[![License: Custom](https://img.shields.io/badge/License-Custom-yellow.svg)](LICENSE)
[![Julia Version](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org/)

**BioToolkit.jl** is a comprehensive, high-performance bioinformatics toolkit for Julia. It provides a unified, Biopython-inspired API for a vast range of workflows: from sequence analysis and structural biology to population genetics, single-cell omics, and GPU-accelerated genomics.

Unlike wrappers around existing tools, BioToolkit features **original, optimized implementations** of core algorithms, delivering significant speedups over traditional Python pipelines and competitive performance with specialized Julia libraries.

## Key Features

*   **High-Performance Core**: Original implementations of parsers and algorithms. Up to **300x faster** than Biopython on core workflows, and **5x faster** than alternative Julia parsing libraries for formats like BED and GFF.
*   **Unified Namespace**: No need to juggle 10 different packages. Sequence I/O, phylogenetics, and structural biology live in one cohesive namespace.
*   **GPU Acceleration**: Native CUDA kernels for k-mer counting, coverage histograms, and massive distance matrices.
*   **Modern Architecture**: Lightweight core with optional extensions for plotting (Makie/Plots) and probabilistic programming (Turing).
*   **Broad Scope**:
    *   **Sequence & Protein Analysis**: GC content, translation, codon usage, melting temps, ProtParam stats.
    *   **Structural Biology**: PDB/mmCIF parsing, RMSD/Kabsch superposition, SASA, contact maps.
    *   **Phylogenetics & PopGen**: Neighbor-joining, consensus trees, $F$-statistics, GWAS scanning.
    *   **Omics**: Differential expression, single-cell analysis, epigenomics, metabolomics, and systems-biology modeling.
*   **DESeq2 + edgeR-style Normalization**: Complete DESeq2 workflow for differential expression analysis, including median-ratio size factors, edgeR-style TMM normalization, dispersion modeling, and Wald testing. See [Examples/deseq2_example.jl](Examples/deseq2_example.jl) for details.
*   **limma-style Empirical Bayes**: Pure-Julia moderated linear modeling with fitFDist/squeezeVar/eBayes-style variance moderation via `limma_fit` in `SystemsBio`.
*   **Single-Cell Interop**: Seurate/Scanpy-like API for H5AD I/O, marker-based cell type annotation, and ambient RNA removal. See [Examples/counts_matrix_projection.jl](Examples/counts_matrix_projection.jl) for details.

## Performance Highlights

BioToolkit is designed for speed. Benchmarks performed on Julia 1.12 vs Python 3.14 (Biopython) and native Julia alternatives.

| Workflow | BioToolkit (ms) | Baseline (ms) | Speedup |
| :--- | :--- | :--- | :--- |
| **Hamming Distance** (19.5 kbp) | 0.79 | Biopython: 259.6 | **328x** |
| **K-mer Frequency** (CUDA) | 0.14 | Biopython: 302.6 | **2161x** |
| **BED Parsing** | 7.7 | BED.jl: 38.2 | **5.0x** |
| **GenBank Parsing** (500 recs) | 3.8 | Biopython: 21.4 | **5.7x** |
| **Phylogenetics** (NJ Tree, 200 taxa) | 11.0 | Biopython: 1848.5 | **168x** |
| **Structural Superposition** | 0.06 | Bio.PDB: 0.75 | **12.6x** |

*See the full [Benchmark Report](benchmark_report.md) for details.*

## Installation

BioToolkit is currently unregistered. You can install it directly from the Git repository:

```julia
using Pkg
# Add the package via URL
Pkg.add(url="https://github.com/Aditya747S/BioToolkit.jl.git")
```

## Quick Start

### Sequence Analysis
```julia
using BioToolkit

seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

# Basic analytics
BioToolkit.gc_content(seq)               # => 0.58
BioToolkit.translate_dna(seq)            # => "MAIVMGR*KGAR*"

# Fast k-mer counting
BioToolkit.kmer_frequency(seq, 3)
```

### Genomic Intervals & I/O
```julia
# High-performance file parsing
records = BioToolkit.read_fasta("genome.fasta")

# Interval algebra (1-based closed intervals)
table = BioToolkit.load_arrow_table("variants.arrow")
subset = BioToolkit.filter_region(table, "chr1", 1_000_000, 1_100_000)

# Coverage & BigWig export
hist = BioToolkit.coverage_histogram(table, "chr1", 1000)
BioToolkit.write_bigwig(hist, "coverage.bw")
```

### Structural Biology
```julia
# Parsing and Geometry
structure = BioToolkit.read_pdb("1a1a.pdb")
coords = BioToolkit.coordinate_matrix(structure)

# Superposition & RMSD (Kabsch algorithm)
ref_struct = BioToolkit.read_pdb("reference.pdb")
result = BioToolkit.superpose(structure, ref_struct)
BioToolkit.rmsd(result)  # => 1.234 (Å)
```

### Population Genetics
```julia
# GenePop workflows
gp = BioToolkit.read_genepop("data.gen")
freqs = BioToolkit.allele_frequencies(gp)
hw_results = BioToolkit.hardy_weinberg_test(gp)
```

## Optional Backends

BioToolkit keeps the core lean. Heavy dependencies for plotting and Bayesian analysis are loaded via Julia extensions:

*   **Plotting**: `import BioToolkitPlotsExt` (Plots.jl) or `import BioToolkitMakieExt` (Makie).
*   **Bayesian**: `import BioToolkitTuringExt` for probabilistic modeling hooks.

Simply installing and importing these packages enables the extra functionality in BioToolkit automatically.

## Documentation

*   **[API Overview](BioToolkit.md)**: Detailed book-level overview of modules.
*   **[Benchmark Report](benchmark_report.md)**: Detailed performance comparisons.
*   **[Module-Specific Docs](src/)**: In-depth documentation for each module (e.g., `structure.md`, `singlecell_interop.md`).

## Project Status

BioToolkit is currently in **Alpha/Active Development**. There can be rough edges. The core API is mostly stable, but some modules are still being fleshed out. Performance optimizations are ongoing, and more workflows will be added over time.

### Development Setup
For development, clone the repo and instantiate dependencies:

```bash
git clone https://github.com/your-username/BioToolkit.jl.git
cd BioToolkit.jl
julia --project
```

```julia
# Inside the Julia REPL
using Pkg
Pkg.instantiate()
Pkg.activate(".")
```

If possible use [Revise.jl](https://timholy.github.io/Revise.jl/stable/) for live code updates during development:

```julia
using Revise
using BioToolkit
```

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue on the repository to discuss or make a proper pull request with added functionality, tests and updated markdown documentation.

## License

BioToolkit is licensed under custom license.
Please see the [LICENSE](LICENSE) file for details.