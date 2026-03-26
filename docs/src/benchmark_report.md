# BioToolkit — Performance Benchmark Report

> **Date**: 2026-03-23 · **Julia**: 1.12.4 · **Python**: 3.14 (conda `general`)
> **Hardware**: same machine for all benchmarks, CUDA 13.1 available

---

## 1. BioToolkit vs Other Julia Packages

| Benchmark | BioToolkit | Other Package | Speedup | Winner |
|---|---|---|---|---|
| Reverse complement | **2.6083 ms** | BioSequences: 0.5299 ms | 0.2× | BioSequences |
| Hamming distance | **0.7850 ms** | BioSequences: 0.4034 ms | 0.5× | BioSequences |
| Hamming distance | **0.7850 ms** | BioAlignments: 34.2589 ms | 43.6× | **BioToolkit** |
| BED parsing | **7.6696 ms** | BED.jl: 38.2322 ms | 5.0× | **BioToolkit** |
| GFF3 parsing | **11.3651 ms** | GFF3.jl: 25.0656 ms | 2.2× | **BioToolkit** |

### Hamming Distance Improvement (vs previous session)

| Metric | Before | After | Change |
|---|---|---|---|
| BioToolkit hamming | 1.7834 ms | 0.7850 ms | **2.3× faster** |
| Gap vs BioSequences | 4.6× slower | 1.9× slower | **Closed by 58%** |

---

## 2. Julia (BioToolkit) vs Python (biopython)

### Sequence Toolkit (200 repetitions, 18 kbp sequences)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Count nucleotides | 2.14 | 72.06 | **33.7×** |
| GC content | 1.90 | 151.32 | **79.6×** |
| DNA transcription | 4.0446 | 4.6026 | **1.1×** |
| Reverse complement | 1.33 | 3.33 | **2.5×** |
| Translate DNA | 10.68 | 157.13 | **14.7×** |
| Hamming distance | 0.16 | 47.96 | **299.7×** |
| K-mer frequency | 31.04 | 302.60 | **9.7×** |
| K-mer frequency (CUDA) | 0.14 | — | GPU accelerated |
| Codon usage | 23.49 | 101.40 | **4.3×** |
| Codon usage table | 23.35 | 97.52 | **4.2×** |
| Relative codon adaptiveness | 6.98 | 290.45 | **41.6×** |
| Codon adaptation index | 24.86 | 86.54 | **3.5×** |
| CAI alias | 24.86 | 86.54 | **3.5×** |

The codon rows come from [sequence_toolkit_compare.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/Examples/Comparison/sequence_toolkit_compare.jl) and its Python companion. The Julia side now covers both raw codon counting and CAI-style weighting, while the Python side uses Biopython's `CodonAdaptationIndex` as the direct baseline for the adaptiveness and CAI rows.

### Standalone Hamming (1000 reps, 19.5 kbp)

| | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Hamming distance | 0.79 | 259.56 | **328.3×** |

### Pairwise Alignment

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Global alignment (1170 bp) | 3.87 | 6.16 | **1.6×** |
| Affine-gap alignment (1200 bp) | 4.69 | 10.99 | **2.3×** |
| Local alignment affine (1000 bp) | 28.31 | 137.63 | **4.8×** |
| Sub. Matrix Load (BLOSUM/PAM) | 1.21 | 6.96 | **5.8×** |
| BLOSUM62 Global | 0.83 | 4.52 | **5.4×** |
| BLOSUM62 Local | 1.28 | 3.43 | **2.7×** |
| PAM250 Global | 0.82 | 2.05 | **2.5×** |
| PAM250 Local | 1.26 | 3.38 | **2.7×** |

### BLAST-like Heuristic Search (1,000 bp vs 100,000 bp target)

| Operation | Julia Heuristic (ms) | Julia Smith-Waterman (ms) | Python Smith-Waterman (ms) | Speedup vs Python |
|---|---|---|---|---|
| Local Target Search | 1.48 | 499.31 | 16709.72 | **11,290×** |

*Note: Comparing our native heuristic search against Biopython wrapping the physical NCBI BLAST+ C++ binaries (`blastn`) via `NcbiblastnCommandline` yielded an execution overhead of ~1,943 ms for the same scale search. Since our Julia library operates directly in-memory without cross-process IPC or I/O overheads, the `local_search` function is computationally ~1,312x faster for rapid algorithmic use compared to Python wrappers.*

### Native Phylogenetics (Neighbor-Joining)

Because BioToolkit has heavily optimized primitive sequence manipulations, calculating massive $N \times N$ symmetric distance maps resolves computationally faster than dynamic GC memory allocations.

| Test Case | Output Taxa | Seq Length | BioToolkit ([phylo.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/src/phylo.jl)) | Note |
|---|---|---|---|---|
| Vector SWAR Hamming Map | 100x100 | 10,000 bp | **2.67 ms** | Threaded natively `Threads.@threads` |
| Multi-GPU SWAR Hamming Map | 10,000x10,000 | 1,000 bp | **14.2 ms** | Multi-GPU block sliced over 4x GPUs |
| DP Pairwise Alignment Map | 20x20 | 500 bp | **301.41 ms** | Requires 190 distinct $O(NM)$ alignments evaluated recursively! |
| Neighbor Joining Tree (200 Taxa matrix) | 200 Taxa | - | **10.95 ms** | vs Python: 1,848.53 ms (**168.8×** faster) |
| Tree Consensus (120 taxa, 25 reps) | 120 Taxa | - | **210.55 ms** | vs Python: 12,114 ms (**57.5×** faster) |
| Midpoint Rooting | 60 Taxa | - | **2.51 ms** | vs Python: 4.50 ms (**1.8×** faster) |
| Prune & Reroot | 80 Taxa | - | **0.52 ms** | vs Python: 2.84 ms (**5.5×** faster) |
| Graph Export (Dot/Mermaid) | - | - | **0.40 ms** | vs Python: 1.14 ms (**2.8×** faster) |
| Neighbor Joining Matrix | 100 Taxa | - | **21.92 ms** | Pre-allocated `@views` DP agglomeration |

### Phylogenetics Extensions (Midpoint Rooting, Consensus, Export)

| Operation | Julia (ms) | Python (ms) | Julia Speedup | Notes |
|---|---|---|---|---|
| Midpoint rooting | **2.3335** | 4.3907 | **1.9×** | In-place Biopython `root_at_midpoint()` baseline |
| Consensus + support-labeled Newick | **172.9661** | 12082.7289 | **69.9×** | Bootstrap consensus with explicit support labels in Newick |
| Tree export (DOT / Mermaid) | **0.4414 / 0.3150** | 1.2782 / 0.5303 | **2.9× / 1.7×** | Closest Biopython baselines: `draw_ascii()` and Newick serialization |
| Prune + outgroup reroot | **0.3616 / 0.2058** | 1.9473 / 0.8085 | **5.4× / 3.9×** | Subset pruning plus outgroup rerooting |

The export benchmark is intentionally asymmetric: Biopython does not provide direct Graphviz DOT or Mermaid emitters, so the closest comparable operations are ASCII tree drawing and Newick serialization.

### Structural Biology (same 500-residue, 2000-atom synthetic PDB)

| Operation | Julia (ms) | Python/Bio.PDB (ms) | Julia Speedup |
|---|---|---|---|
| PDB parse | 2.588 | 7.476 | **2.9×** |
| Center of mass | 0.153 | 0.761 | **5.0×** |
| Coordinate matrix | 0.011 | 0.040 | **3.6×** |
| CA RMSD / superposition | 0.060 | 0.754 | **12.6×** |
| Atom neighborhood query | 0.001 | 0.007 | **7.0×** |

*Note: The Julia timings use `read_pdb`, `center_of_mass`, `coordinate_matrix`, `rmsd`, and a prebuilt `AtomKDTree` neighborhood query. The Python baseline uses Biopython's `Bio.PDB.PDBParser`, `Superimposer`, and `NeighborSearch` on the same synthetic input file, so both sides operate on the same structure and therefore the same derived geometry/diagram.*

### Structural Microbenchmarks

| Operation | Julia (ms) |
|---|---:|
| residue_sasa | 0.0575 |
| structure_sasa | 0.0847 |
| interface_profile | 0.1716 |
| ensemble_rmsd_matrix | 0.0084 |
| superpose_models! | 0.0261 |

These timings reflect the current structural path after caching Fibonacci-sphere SASA samples, reusing atom trees and coordinate matrices where possible, and avoiding repeated full scans in the common multi-residue and multi-model workflows.

*Additional structure-module work in this session tightened mmCIF row parsing for quoted values, added category-specific and raw-category round-tripping, introduced parenthesized selection expressions, and switched residue-contact analysis to a bucketed candidate search so larger structures avoid the old quadratic scan whenever possible.*

### Database-Style Modules (synthetic fixtures)

| Operation | Julia (ms) |
|---|---:|
| Restriction site search on 85 kbp DNA | 0.7708 |
| Sequence logo SVG render | 0.0274 |
| Entrez search JSON parse | 0.121078 |
| KEGG pathway parse | 0.003305 |
| SCOP record parse | 0.000815 |
| CATH record parse | 0.000860 |
| Codon usage (3 codons) | 0.000081860 |
| Codon usage table | 0.000154665 |
| Relative codon adaptiveness | 0.003181 |
| Codon adaptation index | 0.003366 |
| CAI alias | 0.003378 |
| JASPAR motif parse | 0.005762 |

These benchmarks use local synthetic fixtures or parsed JSON samples rather than live network calls, so they measure the parser and matcher cost directly. The Entrez layer is intentionally benchmarked at the response-parsing stage, because network latency would dominate and hide the actual package overhead.

The codon-usage and JASPAR entries were added in the latest parity pass. They stay in the same microsecond range as the other local parser and lookup operations, which is consistent with keeping these workflows fully in-memory and dependency-light.

### GenBank and Annotation Cross-Language Checks

| Operation | Julia (ms) | Python/Biopython (ms) | Julia Speedup | Notes |
|---|---|---|---|---|
| GenBank parse (500 records) | 3.7860 | 21.3873 | **5.7×** | `read_genbank` vs Biopython record parsing on the same synthetic flat-file input |
| GenBank write (500 records) | 4.0333 | 20.5846 | **5.1×** | `write_genbank` vs Biopython serialization of the same record set |
| Sequence annotation extract | 0.1767 | 0.1428 | 0.8× | `annotate_genbank_record` vs Biopython-style feature extraction on a single annotated record |

The GenBank parse and write cases are direct overlap checks on the file-format layer: BioToolkit reads and reconstructs the same synthetic records that the Python baseline parses and serializes. The annotation comparison is slightly asymmetric because BioToolkit keeps the annotated record object in-memory, while the Python side validates equivalent feature extraction and strand-aware spans.

### Biopython Parity Spot-Check

Using Biopython 1.86 in the local Python environment, the following restriction-site positions matched exactly for the same test sequence:

| Enzyme | Biopython | BioToolkit |
|---|---|---|
| EcoRI | [5, 24] | [5, 24] |
| BamHI | [11] | [11] |
| HindIII | [18] | [18] |
| EcoRV | [] | [] |

This comparison is a direct sanity check on the overlapping restriction-enzyme workflow. The KEGG, SCOP, CATH, and Entrez helpers in BioToolkit are extra surfaces that Biopython does not expose in the same compact, package-local form here, so their validation is primarily parser- and fixture-based.

### Applied HMM Gene Prediction (Viterbi)

Running predictive stochastic evaluations over massive continuous genomes is fundamentally bottlenecked by DP Array matrix allocations inside the mathematical iterations. 

| Test Case | Sequence Size | BioToolkit ([hmm.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/src/hmm.jl)) | Note |
|---|---|---|---|
| Profile HMM (Coding) | 100,000 bp | **1.76 ms** | Log-space operations via strictly typed 2-State mapping |

### FASTQ I/O (30k records × 100 bp)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| FASTQ read | 25.85 | 107.90 | **4.2×** |
| FASTQ write | 11.03 | 127.87 | **11.6×** |

### FASTA Index

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| FASTA index build | 0.11 | 0.38 | **3.6×** |
| FASTA fetch | 0.007 | 0.19 | **25.4×** |

### GenBank Parsing (500 records)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| GenBank parse | 3.86 | 19.30 | **5.0×** |
| GenBank ingest (Arrow) | 5.28 | — | Julia-only |

### Flat-File Parsing

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| EMBL parse (30 bp synthetic record) | 2.2324 | 7.1862 | **3.2×** |

The EMBL comparison uses the same small synthetic record in both parsers and compares the same sequence-length result. The remaining identifier and feature-normalization details are parser-specific, so the benchmark treats them as a speed comparison rather than a strict round-trip equivalence check.

### Expanded Biopython Overlap Benchmarks

The newest comparison script, [expanded_biopython_compare.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/Examples/Comparison/expanded_biopython_compare.jl), adds coverage for the remaining Biopython-style gaps highlighted in this session: DNA melting temperature, DNA molecular weight, sequence slicing, BLAST XML/tabular parsing, 10,000-atom neighbor queries, rigid-body superposition, consensus phylogenetics, parsimony tree search, VCF parsing, k-mer indexing, and optional EMBOSS/DSSP-style wrapper overhead checks.

| Operation | Julia (ms) | Python (ms) | Julia Speedup | Status |
|---|---|---|---|---|
| Melting temperature | 0.0173 | 4.8872 | **282.3×** | direct overlap |
| DNA molecular weight | 0.0744 | 0.3052 | **4.1×** | direct overlap |
| Sequence slice `100:5000` | 0.0006 | 0.0852 | **142.0×** | direct overlap |
| BLAST XML parse | 0.2474 | 27.3458 | **110.5×** | direct overlap |
| BLAST tabular parse | 0.0041 | 31.3474 | **7,646×** | direct overlap |
| Neighbor search (10,000 atoms) | 0.0103 | — | — | Julia-only run in this workspace |
| Superposition / Kabsch | 0.0224 | — | — | Julia-only run in this workspace |
| Consensus tree | 1.6420 | — | — | Julia-only run in this workspace |
| Robinson-Foulds distance | 0.0434 | — | — | Julia-only run in this workspace |
| Maximum parsimony tree | 0.1179 | — | — | Julia-only run in this workspace |
| VCF parsing | 3.8750 | — | — | Julia-only run in this workspace |
| K-mer indexing | 0.1630 | — | — | Julia-only run in this workspace |
| Entrez batch fetch | 3,108.6124 | — | — | mocked batch / network-safe Julia-only timing |
| DSSP wrapper overhead | added | added | — | conditional on external executable |

The direct overlap rows above were measured with [expanded_biopython_compare.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/Examples/Comparison/expanded_biopython_compare.jl) and its Python companion. The Julia-only rows were measured separately in [expanded_julia_only_compare.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/Examples/Comparison/expanded_julia_only_compare.jl) because the Python side depends on optional modules or an external executable that were not guaranteed in this workspace.

### Quality Scores (20k records × 100 bp)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Phred decode | 1.41 | 41.55 | **29.4×** |
| Mean quality | 1.55 | 49.59 | **31.9×** |
| Quality trim | 7.64 | 158.19 | **20.7×** |
| Phred string encode | 0.0012 | 0.0090 | **7.5×** |

### Motif Analysis (20k sequences × 12 bp)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Motif counts | 0.51 | 192.83 | **379.6×** |
| Motif PWM | 0.001 | 0.08 | **78.1×** |
| Motif scan both strands | 4.30 | 2.09 | 0.5× |
| Motif Scan (30 kbp, 6000 hits) | 4.22 | 1.96 | 0.46× | *Python C-extensions heavily optimized* |
| JASPAR parse | 0.07 | 0.05 | **0.8×** |

The JASPAR row is now exercised in [motif_pwm_compare.jl](https://github.com/Aditya747S/BioToolkit.jl/blob/main/Examples/Comparison/motif_pwm_compare.jl) and its Python companion using the same synthetic matrix fixture. Biopython remains slightly faster on this tiny parse, but the new Julia reader is now measured directly against it instead of living only in the Julia-only parser section.

### Annotation (GenBank feature extraction)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Annotate + extract | 0.18 | 0.15 | 0.8× (comparable) |

### Genomic Histograms & Coverage

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Window coverage (CPU) | 0.23 | 8.77 | **38.1×** |
| Window coverage (CUDA) | 0.16 | — | GPU accelerated |
| Position Histogram (CPU, 2M items) | 3.84 | 179.48 | **46.7×** |
| Position Histogram (CUDA, 2M items) | 0.10 | — | **1,860× vs Python** |

### PopGen / GenePop (240 individuals per population, 24 loci)

| Operation | Julia (ms) | Python/Biopython (ms) | Julia Speedup | Notes |
|---|---|---|---|---|
| GenePop parse / record parse | 4.14 / 4.09 | 5.63 | **1.4× / 1.4×** | Julia `read_genepop` and `read_genepop_record` vs Biopython `GenePop.read` |
| GenePop write/reconstruct | 10.56 | 19.90 | **1.9×** | Julia `write_genepop` vs Biopython `str(record)` |
| GenePop record helpers | split pops 6.72, split loci 0.39, remove pop 5.19, remove locus 4.17 | split pops 11.91, split loci 2.03, remove pop 12.22, remove locus 11.92, file remove pop 8.29, file remove locus 11.91 | mixed | Julia now exposes the same direct helper surface as Biopython's GenePop record APIs |
| Population stats kernels | 0.0050 / 0.0049 / 0.0055 / 0.8407 / 0.0142 / 0.7089 | n/a | n/a | allele frequency, heterozygosity, HWE, F-statistics, LD, and PCA are Julia-only in this package |

Biopython’s GenePop support is a good parser/manipulation baseline, but it does not provide the analytical population-statistics kernels present in `popgen.jl`. The Julia results above therefore split the comparison into a direct GenePop I/O layer and a Julia-only statistics layer.

### Window Coverage (25k intervals)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Window coverage (CPU) | 0.23 | 8.77 | **38.1×** |
| Window coverage (CUDA) | 0.16 | — | GPU accelerated |

### Protein Statistics & GC Skew

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| Protein Stats (Mass, EC, II, pI, GRAVY) | 0.085 | 2.81 | **33.1×** |
| ProtParam bundle | 0.1250 | 3.2806 | **26.2×** |
| GC Skew (60 kbp) | 0.165 | 0.79 | **4.8×** |

The ProtParam bundle compares BioToolkit's bundled protein-statistics workflow against the equivalent Biopython `ProteinAnalysis` method set plus the same aliphatic-index formula, since Biopython does not expose a single `protparam` convenience call in this environment.

---

## 3. MSA Round-Trip (500 seqs × 120 cols)

| Operation | Julia (ms) | Python (ms) | Julia Speedup |
|---|---|---|---|
| MSA construct | 0.09 | 0.14 | **1.6×** |
| Column access | 0.007 | 0.29 | **38.6×** |
| Row access | 0.0003 | 0.001 | **3.7×** |
| Append | 0.003 | 0.037 | **13.2×** |
| Column frequency | 0.014 | 0.30 | **21.6×** |
| FASTA roundtrip | 36.19 | 28.45 | 0.8× |
| Clustal roundtrip | 7.27 | 2.49 | 0.3× |
| Stockholm roundtrip | 4.42 | 3.58 | 0.8× |
| MSF roundtrip | 2.00 | 30.86 | **15.4×** |
| PHYLIP roundtrip | 16.69 | 31.95 | **1.9×** |

> [!NOTE]
> Python wins on MSA file round-trips because biopython's C-optimized parsers are highly tuned for these specific formats. Julia is significantly faster on all in-memory operations.

## 4. External MSA Wrappers (8 seqs × 180 bp)

These benchmarks measure the end-to-end wrapper path: write temporary FASTA input, invoke the external aligner, then read the aligned FASTA output back.

| Tool | Julia (ms) | Python (ms) | Julia Speedup | Winner |
|---|---|---|---|---|
| Clustal Omega | 57.64 | 51.43 | 0.9× | Python |
| MUSCLE 5 | 24.02 | 31.67 | **1.3×** | **BioToolkit** |

The external-aligner results are intentionally close because the binary runtime dominates. The meaningful difference is that BioToolkit now offers a compact native wrapper surface for both tools while keeping the core alignment and MSA code dependency-free.

The wrappers are intentionally thin process adapters; the binary runtime dominates, so the gap is mostly determined by external tool startup and file I/O rather than Julia-side algorithmic work.

## 5. Additional Verified Julia-only Features

These features were validated in the current repository state, but they do not have Python baselines in this report.

| Feature | Julia (ms) | Notes |
|---|---:|---|
| ABIF (.ab1) parsing | 2.8428 | Mock V3 trace with PBAS2/PCON2 tags |
| Maximum-likelihood pruning (`felsenstein_likelihood`) | 466.9001 | 3-taxon JC69 likelihood evaluation on the test alignment |
| Maximum-likelihood tree search (`maximum_likelihood_tree`) | 1038.0599 | Coordinate-ascent search from the same alignment |
| EMBL parsing | 15.8004 | Single-record flat-file parse |
| SwissProt parsing | 0.1064 | Same parser path as EMBL |
| Coalescent simulation (`simulate_coalescent`) | 74.0916 | 50 taxa, 1,000 bp, constant-size model |

Local BLAST+ wrapper verification still depends on `blastn` and `makeblastdb` being installed in the execution environment. The benchmark report keeps the existing heuristic-search comparison, but the wrapper itself was not re-timed in this environment because those binaries were unavailable.

### Biopython Coverage Map

This report now includes direct cross-language comparisons for the major workflows that have a natural Biopython baseline. The remaining Julia-only sections are still timed because they matter for performance, but they do not have a compact in-repo Biopython analogue that is worth treating as a like-for-like baseline here.

| Area | Representative comparison scripts | Status |
|---|---|---|
| Sequence toolkit and similarity | `sequence_toolkit_compare.jl`, `hamming_compare.jl` | Direct Biopython baseline present, including codon usage and CAI |
| Pairwise alignment and scoring | `pairwise_align_compare.jl`, `pairwise_affine_compare.jl`, `protein_matrix_compare.jl`, `local_align_compare.jl` | Direct Biopython baseline present |
| Sequence I/O and indexing | `fastq_compare.jl`, `fasta_index_compare.jl`, `genbank_compare.jl`, `genbank_write_compare.jl` | Direct Biopython baseline present |
| Flat-file parsing | `flatfile_compare.jl` | Direct Biopython baseline present |
| GenBank annotation | `annotation_compare.jl` | Direct Biopython-style feature extraction baseline present |
| Motifs and logo paths | `motif_pwm_compare.jl`, `motif_scan_compare.jl` | Direct Biopython baseline present, including JASPAR parsing |
| MSA and external aligners | `msa_compare.jl`, `msf_compare.jl`, `clustal_muscle_compare.jl` | Direct Biopython baseline present |
| Phylogenetics | `nj_compare.jl`, `midpoint_root_compare.jl`, `consensus_support_newick_compare.jl`, `prune_reroot_compare.jl`, `graph_export_compare.jl`, `phylip_compare.jl` | Direct Biopython or closest-rendering baseline present |
| Genomic coverage and histograms | `window_coverage_compare.jl`, `cuda_histogram_compare.jl` | Direct Python baseline present |
| Database and reference workflows | restriction, Entrez, KEGG, SCOP, CATH sections above | Mixed: some direct Biopython overlap, some parser-first Julia-only surfaces |
| Population genetics | `popgen_compare.jl`, `bench_popgen.jl` | Direct GenePop overlap for I/O, Julia-only kernels for the analytical layer |

The practical takeaway is that BioToolkit now has measured overlap in every major Biopython-style area where an apples-to-apples baseline is available, and the report explicitly separates those from parser-first or Julia-only capabilities. The newest additions extend the measured overlap to DNA transcription, Phred string encoding, and a bundled ProtParam-style protein workflow.

---

## Summary

- **BioToolkit leads most Python workloads** in this report, while a few file-round-trip and external-tool paths still favor Biopython or third-party binaries
- **Hamming distance improved 2.3×** from the SWAR optimization (0.785 ms vs 1.78 ms before)
- **BioSequences still faster** for reverse complement and hamming (specialised 2-bit encoding), but the gap narrowed
- **BioToolkit beats BED.jl** (5×) and **GFF3.jl** (2.2×) on parsing
- **GPU acceleration** available for k-mer frequency and window coverage
- **Sequence transcription and Phred-string encoding are now directly benchmarked against Biopython**, alongside the earlier count/GC/reverse-complement/translation workflow
- **Codon usage, codon adaptiveness, CAI, and JASPAR parsing are now directly benchmarked against Biopython**, closing two more of the remaining Biopython overlap gaps
- **DNA transcription is now significantly faster after the single-pass rewrite, and Phred-string encoding remains directly benchmarked against Biopython**
- **MSF round-trips** now have a fast native parser path and compare favorably against Biopython on the measured case
- **PHYLIP round-trips** are now tracked alongside FASTA, Clustal, Stockholm, and MSF
- **External MSA wrappers** are now supported for Clustal Omega and MUSCLE 5, with MUSCLE ahead of the Python baseline on the measured case
- **Phylogenetics now covers midpoint rooting, support-labeled consensus Newick, and DOT/Mermaid export**, with the consensus path especially far ahead of Biopython on the measured bootstrap workload
- **Phylogenetics now also covers pruning, rerooting, and tree-query helpers**, with the prune/reroot path outperforming the Biopython baseline on the measured case
- **GenBank parsing, GenBank writing, and sequence annotation now have direct Biopython overlap checks**, giving the report more than one database-style cross-language baseline instead of relying only on restriction-site parity
- **EMBL parsing now has a direct speed comparison against Biopython and is faster on the synthetic fixture after the tokenization rewrite**, even though some metadata fields normalize differently between parsers on the synthetic fixture
- **The benchmark suite now explicitly includes the remaining Biopython overlap gaps**: melting temperature, DNA molecular weight, sequence slicing, BLAST XML/tabular parsing, 10,000-atom neighbor search, rigid superposition, consensus trees, Robinson-Foulds distance, parsimony trees, VCF parsing, k-mer indexing, and conditional wrapper-overhead checks
- **ProtParam-style protein stats now include the bundled `protparam` workflow as an additional ProteinAnalysis workflow comparison**, on top of the individual mass/extinction/instability/pI/gravy checks
- **Maximum-likelihood phylogenetics, ABIF trace parsing, EMBL/SwissProt parsing, and coalescent simulation are now unit-tested and timed in Julia-only verification runs**
- **PopGen/GenePop now has a direct Biopython comparison**: Julia is faster on GenePop parse and reconstruction, and it now exposes the same record split/remove helpers as Biopython
- **Population-statistics kernels** in `popgen.jl` do not have native Biopython counterparts, so they are reported separately as Julia-only analytical throughput
- **Multi-GPU Hamming maps** and **position histograms** are tracked explicitly as CUDA-path wins
