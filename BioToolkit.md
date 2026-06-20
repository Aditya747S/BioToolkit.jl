# BioToolkit.jl

**High-Performance, Pure-Julia Computational Biology & Bioinformatics**

BioToolkit.jl is a comprehensive, end-to-end framework for modern bioinformatics written entirely in Julia. By leveraging Julia's multiple dispatch, parametric type system, and just-in-time compilation, BioToolkit eliminates the runtime overhead and type-safety pitfalls of traditional Python/R workflows while matching or exceeding the performance of C/C++ backends.

---

## Core Philosophy

| Principle | Implementation |
|-----------|----------------|
| **Strict Type Safety** | Sequences are parametrically typed by alphabet (`DNASeq`, `RNASeq`, `AASeq`), preventing silent cross-alphabet corruption (e.g., translating a DNA sequence as a protein) at compile time. |
| **Zero-Copy & Cache-Optimized** | Inner loops operate on raw `Vector{UInt8}` using SWAR (SIMD Within A Register) algorithms, 256-byte lookup tables, and prefix-sums to maximize L1 CPU cache hits. |
| **Sparse-First Math** | Single-cell and genomics matrices default to `SparseMatrixCSC`, allowing analyses on million-cell datasets using megabytes instead of gigabytes of RAM. |
| **GPU-Ready** | Transparent CUDA offloading via lazy `@eval` dispatch—pass a `CuArray` to trigger GPU kernels, or pass a standard `Array` for optimized multi-threaded CPU BLAS. Zero compile-time dependencies on CUDA.jl. |
| **Data-Plot Separation** | Plotting functions return strongly-typed layout structs (`VolcanoPlotResult`, `GeneRenderPlan`), allowing headless computation, caching, or rendering via any backend (Makie, Plots, Web). |

---

## Architecture Overview

BioToolkit is organized into modular, loosely-coupled domains. You can use it as a lightweight sequence parser or as a full multi-omics analysis engine.

```text
┌─────────────────────────────────────────────────────────────┐
│                    VISUALIZATION & BROWSERS                 │
│         BioPlotting.jl │ GenomeBrowser.jl │ HMM.jl          │
├─────────────────────────────────────────────────────────────┤
│                     MULTI-OMICS & CLINICAL                   │
│  DifferentialExpression │ Enrichment │ GWAS │ Clinical.jl  │
├──────────────┬──────────┬──────────────┬────────────────────┤
│   EPIGENOMICS │  SINGLE  │  PHYLOGENY  │  SYSTEMS BIO      │
│ Epigenetics.jl│  CELL    │  Phylo.jl    │  systemsbio.jl    │
│              │  Single  │  PopGen.jl   │  Microbiome.jl    │
│              │  Cell.jl │              │  Metabolomics.jl  │
├──────────────┴──────────┴──────────────┴────────────────────┤
│                   SEQUENCES & STRUCTURE                     │
│  Sequence.jl │ Motif.jl │ Protein.jl │ structure.jl       │
├─────────────────────────────────────────────────────────────┤
│                  ALIGNMENT & SEARCH                         │
│  align.jl │ MSA.jl │ Search.jl │ Proteomics.jl             │
├─────────────────────────────────────────────────────────────┤
│              GENOMIC GEOMETRY & INTERVALS                   │
│  GenomicRanges.jl │ annotation.jl │ GenomicIntervals.jl   │
├─────────────────────────────────────────────────────────────┤
│                      I/O & DATABASES                        │
│  IO.jl (FASTA/Q, VCF, GFF, BAM, GenBank) │ Databases.jl    │
├─────────────────────────────────────────────────────────────┤
│                        FOUNDATIONS                          │
│  biotypes.jl (BioSequence, Rle, IntervalTree, SummarizedExp)│
└─────────────────────────────────────────────────────────────┘
```

---

## Feature Domains

### 1. Foundations (`biotypes.jl`)
The bedrock of BioToolkit. Defines the `BioAlphabet` singleton types (`DNAAlphabet`, `RNAAlphabet`, `AminoAcidAlphabet`) and the core `BioSequence{A}` container. Also provides memory-efficient `Rle` (Run-Length Encoding), `SummarizedExperiment` (coordinated assay + metadata), and an AVL-augmented `IntervalTree` for $O(\log n + k)$ spatial queries.

### 2. Sequence Operations (`Sequence.jl`, `QualityControl.jl`)
Extreme-performance nucleotide manipulation.
* **Translation:** 6-bit codon lookups (`translate_dna`).
* **K-mers:** Base-4 packed `UInt64` hashing for $O(1)$ dictionary lookups (`kmer_frequency`).
* **Distances:** SWAR-optimized Hamming distance comparing 8 bytes per CPU cycle.
* **QC:** $O(N)$ prefix-sum sliding window trimming (`trim_low_quality`), seeded adapter matching (`adapter_trim`), and end-to-end short-circuiting pipelines (`process_sequencing_record`).
* **Genome Access:** Zero-copy random access FASTA via `Mmap.mmap` and `.fai` indexing.

### 3. Alignment & Search (`align.jl`, `MSA.jl`, `Search.jl`)
* **Pairwise:** Needleman-Wunsch, Smith-Waterman, and Gotoh (affine gaps) with space-optimized $O(N)$ trace matrices.
* **Codon Alignment:** Triplet-strict alignment with packed `UInt8` tokens (9x memory reduction).
* **MSA:** Array-like 2D slicing over standard formats (FASTA, CLUSTAL, Stockholm, NEXUS, etc.), consensus generation, and wrappers for CLUSTALO/MUSCLE.
* **BLAST:** Native seed-and-extend using `UInt64` k-mer indexes, bounded X-drop extensions, and typed XML/TSV parsing of NCBI outputs.

### 4. Genomic Geometry (`GenomicRanges.jl`, `annotation.jl`, `GenomicIntervals.jl`)
* **Intervals:** Strict 1-based closed / 0-based half-open coordinate safety. Sweep-line interval subtraction ($O(N \log N)$), intersection, union, and complementation.
* **Annotations:** Full GenBank flat-file parsing, complex location math (`join`, `complement`), and intelligent record slicing that remaps feature coordinates into new frames (including automatic reverse-complement handling).
* **Coverage:** Event-driven sweep-line depth calculation (`SparseCoverageVector`).

### 5. Next-Gen Sequencing I/O (`IO.jl`, `bam.jl`)
* **Parsing:** Byte-level `findnext` scanning for GFF/BED (zero intermediate allocations), robust VCF quality mapping (`.` $\to$ `missing`), and multi-line GenBank qualifier handling.
* **BAM:** Pure-Julia BGZF decoding, nibble unpacking, CIGAR parsing, auxiliary tag I/O, and standard BAI index generation for $O(\log n)$ random access.
* **Big Data:** Chunked Arrow conversions (`ingest_vcf`, `ingest_gff`) for zero-OOM processing of terabyte-scale datasets.

### 6. Single-Cell & Spatial (`SingleCell.jl`)
A complete Scanpy/Seurat replacement in pure Julia.
* **Preprocessing:** Sparse-triplet library size normalization, Seurat v3 VST feature selection, native SCTransform.
* **Graphs:** Native Leiden and Louvain clustering on Shared Nearest Neighbor (SNN) graphs.
* **Trajectory:** Minimum Spanning Tree (MST) pseudotime, and scVelo-style dynamical RNA velocity using `DifferentialEquations.jl`.
* **Spatial:** Moran's I spatial autocorrelation and spatially variable gene detection.
* **Multi-Modal:** Weighted Nearest Neighbors (WNN) integration (RNA + ATAC).
* **Interoperability:** Full read/write parity with Scanpy `.h5ad` (AnnData) format.

### 7. Epigenomics (`Epigenetics.jl`)
* **Peak Calling:** Event-driven sparse coverage, local Poisson background estimation (MACS2-style), and BH q-value assignment.
* **Methylation:** Cytosine binning, Beta-binomial grid-search dispersion estimation, and Likelihood Ratio Tests (LRT) for DMLs.
* **3D Genomics:** Directionality Index (DI) calculation and Hidden Markov Model (HMM) TAD detection.
* **Footprinting:** Flank vs. center signal depletion for TF binding.

### 8. Multi-Omics Statistics (`DifferentialExpression.jl`, `Enrichment.jl`, `GWAS.jl`, `Clinical.jl`)
* **Differential Expression:** A pure-Julia, zero-dependency reimplementation of DESeq2 (Median-of-ratios, Cox-Reid penalized dispersion, IRLS GLM, Wald/LRT, apeglm LFC shrinkage).
* **Enrichment:** Fisher's Exact Test with Haldane-Anscombe OR correction, Benjamini-Hochberg FDR, and hierarchical redundancy pruning (Elim algorithm).
* **GWAS:** PLINK binary parsing, QR-projected linear scans, spectral LMM transformation, LD clumping, LDpred ridge PRS, and DerSimonian-Laird meta-analysis. Transparent CUDA offloading.
* **Clinical:** Kaplan-Meier, Log-Rank, Cox Proportional Hazards (formula interface), Time-Dependent ROC, Competing Risks (CIF), Neural Cox, and TCGA GDC API ingestion.

### 9. Evolution & Phylogenetics (`Phylo.jl`, `PopGen.jl`)
* **Phylogenetics:** Newick, PhyloXML, NEXUS, and NeXML I/O. NJ & UPGMA tree construction. Maximum Likelihood (JC69, K80, HKY85) via Felsenstein's Pruning. Robinson-Foulds distance, bootstrapping, and consensus methods.
* **Population Genetics:** Strict parametric types (`Locus{T}`). Weir & Cockerham F-stats, AMOVA, Hardy-Weinberg exact tests, LD decay, EHH/iHS/XP-EHH, ABBA-BABA (D/F-statistics), Wright-Fisher/Coalescent simulations.

### 10. Structure & Proteomics (`structure.jl`, `Proteomics.jl`, `Protein.jl`, `Motif.jl`)
* **Structure:** PDB/mmCIF parsing, Kabsch superposition, RMSD, interface/residue contacts, B-factor flexibility, DSSP integration, and Ramachandran analysis.
* **Proteomics:** Wavelet peak picking, Dynamic Time Warping (ObiWarp), QRILC missing value imputation, sparse PLS-DA, and graph-based de novo sequencing.
* **Proteins:** 256-byte table-driven ProtParam calculations (MW, pI, GRAVY, Instability Index) in a single pass.
* **Motifs:** Position Weight Matrix (PWM) math, $O(L)$ 256-byte lookup scanning (dual-strand), de novo seed-and-extend discovery, and pure-SVG sequence logo generation.

### 11. Systems Biology (`systemsbio.jl`, `Microbiome.jl`, `Metabolomics.jl`)
* **Networks:** WGCNA-style soft-thresholding, topological overlap, module detection, and module eigengenes.
* **Microbiome:** Compositional transforms (CLR/ILR), alpha/beta diversity (Faith's PD, Unifrac), ANCOM, Songbird multinomial regression, and co-occurrence networks.
* **Metabolomics:** Bayesian source tracking (Turing.jl), async LC-MS streaming, and delegation to the Proteomics DA engine.

### 12. Visualization & Browsers (`BioPlotting.jl`, `GenomeBrowser.jl`)
* **Plots:** Data-plot separated Volcano, MA, Manhattan, QQ, Forest, and Clustered Heatmap generation.
* **Browsers:** Backend-agnostic layout math engine computing exact track coordinates, row packing, and coverage bins without drawing anything (renderable via Makie, Plotly, etc.).

### 13. Infrastructure (`HMM.jl`, `CUDA.jl`, `Databases.jl`)
* **HMM:** Log-space Viterbi, Forward, and Backward algorithms with $O(1)$ emission lookups.
* **GPU:** Lazy-loaded CUDA kernels for Hamming distances, base-5 k-mer hashing, and multi-GPU phylogenetic distance matrices.
* **Databases:** NCBI E-utilities, MEDLINE/PubMed parsing, EMBOSS CLI wrappers, KEGG flat-file parsing, and Mermaid.js pathway graph export.

---

## Quick Start

### Installation
```julia
using Pkg
Pkg.add(url="https://github.com/Aditya747S/BioToolkit.jl.git")
```

### End-to-End Workflow Example
```julia
using BioToolkit
using BioToolkit.DifferentialExpression
using BioToolkit.BioPlotting

# 1. Load raw counts (or read from a file)
counts = CountMatrix(my_int_matrix, gene_ids, sample_ids)
coldata = DataFrame(condition=["Ctrl", "Ctrl", "Treat", "Treat"])

# 2. Run the DESeq2 pipeline
dds = DESeqDataSet(counts, coldata, ~condition)
dds = DESeq(dds; test=:Wald)
res_df = results(dds)

# 3. Shrink log2 fold changes
shrunk_res = shrink_lfc(res_df)

# 4. Plot the results
vol_plot = volcano_plot(shrunk_res; title="Treatment vs Control", save_path="volcano.png")
```

---

## Migrating from Python / R

BioToolkit is designed as a drop-in replacement for fragmented ecosystems like `Biopython`, `Scanpy`, `Seurat`, and `Bioconductor`. 

| Task | Python / R | BioToolkit.jl |
|------|------------|---------------|
| **Typed Sequences** | `Seq("ATCG")` (String) | `DNASeq("ATCG")` (Parametric Type) |
| **Alignments** | `Bio.Align.substitution_matrices.load("BLOSUM62")` | `named_substitution_matrix("BLOSUM62")` |
| **Genomic Intervals** | `pybedtools.BedTool("a.bed")` | `build_collection(read_bed("a.bed"))` |
| **Differential Exp.** | `scanpy.tl.rank_genes_groups` | `DESeq(dds)` (Pure DESeq2 impl.) |
| **Single-Cell** | `scanpy.pp.normalize_total` | `normalize_counts(exp)` |
| **GWAS** | `statsmodels.GLM` | `gwas_linear_scan(geno, pheno)` |
| **Phylogenetics** | `ete3.Tree("newick")` | `parse_newick("newick")` |

*Unlike Python/R, BioToolkit uses strict compiler-level type enforcement. Passing an `AASeq` to a DNA aligner results in a compile-time `MethodError`, eliminating a massive class of silent data-corruption bugs.*

---

## Detailed API Documentation

BioToolkit is extensively documented. For specific type signatures, keyword arguments, and algorithmic details, please refer to the module-specific API references:

*   [`biotypes.jl`](docs/biotypes.md) — Alphabets, BioSequence, Rle, SummarizedExperiment
*   [`align.jl`](docs/align.md) — Pairwise & Codon Alignment, Scoring Matrices
*   [`longread.jl`](docs/longread.md) — Long-read Seeding, Banded Alignment, SV Calling
*   [`annotation.jl`](docs/annotation.md) — Genomic Locations, Feature Slicing, Variant Annotation
*   [`bam.jl`](docs/bam.md) — Pure-Julia BAM I/O & BAI Indexing
*   [`BioPlotting.jl`](docs/bioplotting.md) — Volcano, Heatmaps, Manhattan, QQ Plots
*   [`Clinical.jl`](docs/clinical.md) — Survival Analysis, MAF, TCGA Integration
*   [`coevolution.jl`](docs/coevolution.md) — DCA Contacts, Reweighting, Contact-guided Folding
*   [`Databases.jl`](docs/databases.md) — Entrez, Restriction Enzymes, EMBOSS, KEGG
*   [`DifferentialExpression.jl`](docs/de.md) — DESeq2 Pipeline, VST, Batch Correction
*   [`Enrichment.jl`](docs/enrichment.md) — ORA, Fisher's Exact, Elim Algorithm
*   [`Epigenetics.jl`](docs/epigenetics.md) — Peak Calling, Methylation, Hi-C, scATAC
*   [`GenePrediction.jl`](docs/geneprediction.md) — HMM Ab Initio Gene Finding
*   [`GenomicRanges.jl`](docs/genomicranges.md) — Interval Algebra, Overlaps, Set Operations
*   [`GWAS.jl`](docs/gwas.md) — Linear/Mixed Models, PRS, Meta-Analysis
*   [`HMM.jl`](docs/hmm.md) — Viterbi, Forward, Backward Algorithms
*   [`IO.jl`](docs/io.md) — FASTA/Q, VCF, GFF, GenBank, Arrow Streaming
*   [`Metabolomics.jl`](docs/metabolomics.md) — Source Tracking, DA, Streaming LC-MS
*   [`Microbiome.jl`](docs/microbiome.md) — ANCOM, Songbird, UniFrac, Co-occurrence
*   [`Motif.jl`](docs/motif.md) — PWM, Scanning, De Novo Discovery, Logos
*   [`MSA.jl`](docs/msa.md) — Multiple Sequence Alignment I/O & Slicing
*   [`pipeline.jl`](docs/pipeline.md) — DAG Workflows, Caching, Template Nodes
*   [`Phylo.jl`](docs/phylo.md) — Tree I/O, ML, Bootstrapping, Consensus
*   [`PopGen.jl`](docs/popgen.md) — F-Stats, HWE, LD, Coalescent, Admixture
*   [`Protein.jl`](docs/protein.md) — ProtParam (MW, pI, GRAVY)
*   [`Proteomics.jl`](docs/proteomics.md) — MS/MS, Imputation, sPLS-DA, De Novo
*   [`QualityControl.jl`](docs/qc.md) — Trimming, Adapter Removal, FastQC stats
*   [`Search.jl`](docs/search.md) — Native BLAST, Seed-and-Extend
*   [`Sequence.jl`](docs/sequence.md) — Translation, K-mers, PCR, CRISPR
*   [`SingleCell.jl`](docs/singlecell.md) — Leiden, Velocity, Spatial, H5AD
*   [`spatial.jl`](docs/spatial.md) — Reference Profiles, RCTD-style and Cell2Location-style Deconvolution
*   [`structure.jl`](docs/structure.md) — PDB, DSSP, Superposition, Interfaces
*   [`systemsbio.jl`](docs/systemsbio.md) — WGCNA, Limma, Network Inference

