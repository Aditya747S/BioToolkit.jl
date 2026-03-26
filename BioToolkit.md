# BioToolkit — Documentation Overview
[BioToolkit](src/BioToolkitBody.jl) is a Julia bioinformatics toolkit with a broad, Biopython-inspired public API. The package combines sequence processing, file I/O, annotations, search, structural biology, phylogenetics, population genetics, and database/reference workflows in one namespace.

This document is the concise book-level overview for the detailed module notes under [BioToolkit.jl/src](BioToolkit.jl/src) and [ext](ext). It points to the source index and optional backend chapters rather than repeating every module-level detail.

Source map:
- [package_overview.md](package_overview.md) is the broader package summary that groups the public API by workflow.
- [ext/BioToolkitPlotsExt.md](ext/BioToolkitPlotsExt.md), [ext/BioToolkitMakieExt.md](ext/BioToolkitMakieExt.md), [ext/BioToolkitCairoMakieExt.md](ext/BioToolkitCairoMakieExt.md), and [ext/BioToolkitTuringExt.md](ext/BioToolkitTuringExt.md) document the optional backend and Bayesian integration layers.

## 1. Core Record and Feature Types
BioToolkit includes lightweight record and feature containers for genomic, trace, and structural workflows.

* Variant records: `VariantEvent`, `VariantTextRecord`
* Interval and annotation records: `BedRecord`, `GffRecord`, `GenBankFeature`, `GenBankRecord`, `GenBankArrowRecord`, `FastqRecord`, `SeqRecordLite`
* Feature-location types: `AbstractFeatureLocation`, `FeatureLocationLite`, `CompoundFeatureLocation`, `SeqFeatureLite`, `AnnotatedSeqRecord`
* Feature utilities: `feature_spans`, `feature_bounds`, `feature_slice`, `feature_sequence`, `feature_extract`, `select_features`, `features_at`, `features_overlapping`, `feature_summary`
* Trace support: `SangerTrace`

## 2. File I/O, Parsing, and Writing
The I/O layer covers the main sequence, trace, flat-file, annotation, and table formats used throughout the package.

* Sequence files: `read_fasta`, `read_fastq`, `write_fastq`
* Trace and flat files: `read_abif`, `read_embl`, `read_swissprot`, `read_blast_xml`, `read_blast_tabular`, `read_medline`
* Genomic formats: `read_vcf`, `write_vcf`, `read_bed`, `write_bed`, `read_gff`, `write_gff`, `read_genbank`, `write_genbank`, `annotate_genbank_record`, `annotate_genbank_records`
* Binary genomics: `BamReference`, `BamHeader`, `BamCigarOp`, `BamRecord`, `BamFile`, `read_bam`, `write_bam`, `write_bigwig`
* Arrow and schema support: `arrow_schema`, `load_arrow_table`, `write_arrow_table`, `encode_chromosome`, `encode_base`, `encode_identifier`, `compact_variant_event`, `parse_bed_record`, `parse_gff_record`, `parse_genbank_record`, `ingest_vcf`, `ingest_bed`, `ingest_gff`
* FASTA indexing: `fasta_index`, `fetch_fasta_sequence`, `FastaIndexRecord`
* NCBI and local search entry points: `run_blast`, `qblast`, `parse_blast_xml`, `parse_blast_tabular`

BED-to-interval conversion adds 1 to the start coordinate, while GFF already maps directly to the package's internal 1-based closed interval convention.

## 3. Sequence Toolkit
The core sequence API handles validation, transformation, composition, comparison, and utility math.

* Validation and composition: `validate_dna`, `count_nucleotides`, `gc_content`, `gc_skew`, `minimum_skew`, `kmer_frequency`, `bin_positions`
* Transformation: `transcribe_dna`, `melting_temp`, `dna_molecular_weight`, `reverse_complement`, `translate_dna`, `find_orfs`, `protein_search`
* Codon usage and CAI: `codon_usage`, `codon_usage_table`, `relative_codon_adaptiveness`, `codon_adaptation_index`, `cai`
* Comparison and visualization: `hamming_distance`, `dotmatrix`

## 4. Protein Statistics
BioToolkit includes ProtParam-style protein property calculations.

* `protein_mass`
* `extinction_coefficient`
* `instability_index`
* `isoelectric_point`
* `gravy`
* `aliphatic_index`
* `protparam`

## 5. Quality and FASTQ Processing
FASTQ helpers cover PHRED handling, trimming, and filtering.

* `phred_score`, `phred_scores`, `phred_string`
* `mean_quality`
* `trim_low_quality`
* `quality_filter`
* `adapter_trim`
* `sequencing_pipeline`

## 6. Pairwise Alignment and Scoring
Pairwise comparison includes general, codon-aware, and matrix-based alignment helpers.

* Alignment result type: `PairwiseAlignmentResult`
* Generic aligners: `pairwise_align`, `needleman_wunsch`, `smith_waterman`, `local_align`
* Codon aligners: `pairwise_align_codons`, `needleman_wunsch_codons`, `smith_waterman_codons`, `local_align_codons`
* Scoring matrices: `SubstitutionMatrix`, `substitution_matrix`, `named_substitution_matrix`, `available_named_substitution_matrices`, `CodonSubstitutionMatrix`, `codon_substitution_matrix`, `named_codon_substitution_matrix`

## 7. Multiple Sequence Alignment and Motifs
The MSA and motif surface covers alignment containers, column statistics, consensus generation, motif discovery, and logo output.

* Alignment containers: `AbstractMultipleSequenceAlignment`, `MultipleSequenceAlignment`, `multiple_sequence_alignment`, `progressive_multiple_sequence_alignment`
* Alignment utilities: `get_alignment_length`, `alignment_column_counts`, `alignment_column_frequencies`, `alignment_column_symbol`, `alignment_column`, `alignment_symbol_line`, `emboss_consensus`, `consensus_sequence`, `clustal_consensus`, `format_alignment`, `write_alignment`, `read_alignment`
* Motif containers: `MotifCounts`, `MotifFrequencyMatrix`, `MotifPWM`, `MotifHit`, `MotifSite`, `MotifDiscoveryResult`, `MotifOccurrence`, `MotifProfile`
* Motif utilities: `motif_counts`, `motif_frequency_matrix`, `motif_entropy`, `motif_relative_entropy`, `motif_scan`, `motif_scan_both_strands`, `motif_information_content`, `discover_motifs`, `motif_pwm`, `motif_consensus`
* Motif parsers and logos: `read_meme`, `read_alignace`, `read_jaspar`, `sequence_logo`, `sequence_logo_svg`, `motif_logo_svg`
* External MSA wrappers: `clustal_msa`, `muscle_msa`

## 8. Search, Indexing, and BLAST-Like Workflows
Search helpers support local index-based scans and structured BLAST parsing.

* `KmerIndex`
* `build_index`
* `HighScoringPair`
* `BlastXMLRecord`, `BlastXMLHit`, `BlastXMLHSP`
* `BlastTabularRecord`, `BlastTabularHit`, `BlastTabularHSP`
* `local_search`
* `parse_blast_xml`, `read_blast_xml`, `parse_blast_tabular`, `read_blast_tabular`

## 9. Hidden Markov Models
The HMM module provides core inference routines and batch scanning.

* `HMM`
* `viterbi`
* `forward`
* `backward`
* `batch_viterbi!`

## 10. Genomic Queries, Coverage, and Interval Algebra
These utilities work over tabular and Arrow-backed genomic data, and the `GenomicRanges` module adds sorted interval collections for direct overlap and neighborhood queries.

* Tabular/Arrow helpers: `filter_region`, `coverage_histogram`, `window_coverage`, `write_bigwig`
* Interval types: `GenomicRanges`, `GenomicInterval`, `IntervalCollection`, `CoverageSegment`, `SeqInfo`
* Interval builders and readers: `build_collection`, `read_intervals`
* Interval queries: `overlap`, `find_overlaps`, `find_overlaps_parallel`, `nearest`, `find_nearest`, `follow`, `precede`
* Interval transforms: `shift`, `flank`, `resize`, `promoters`, `narrow`, `trim`, `gaps`, `disjoin`
* Interval algebra and coverage: `GenomicRanges.setdiff`, `GenomicRanges.reduce`, `intersect`, `union`, `coverage`

`write_bigwig` accepts either a single coverage vector or a chromosome-to-vector dictionary, so multi-chromosome tracks can be emitted without stitching files together manually.

`DataFrame(intervals)` now produces lowercase `chrom`, `start`, `stop`, and `strand` columns so the output round-trips cleanly through `read_intervals`.

## 11. Differential Expression
The `DifferentialExpression` module provides a native RNA-Seq statistical pipeline for sparse count matrices.

* Core containers: `CountMatrix`, `DEResult`, `GLMSolver`
* Normalization and dispersion: `calc_norm_factors`, `estimate_dispersions`
* Statistical testing: `differential_expression`, `benjamini_hochberg`, `fit_gene_fast!`
* Effect-size and preprocessing helpers: `filter_low_counts`, `shrink_lfc`, `vst`
* Batch correction: `remove_batch_effect`, `combat_correction`, `estimate_surrogates`

The implementation is tuned for synthetic and real count data workflows, using sparse input storage, TMM-style normalization, shrinkage dispersion estimates, threaded per-gene fitting, LFC shrinkage, variance-stabilizing transforms, direct linear batch projection, ComBat-style empirical-Bayes correction, and surrogate-variable estimation.

The batch-correction helpers are exported through the root loader so they are available as part of the public BioToolkit API.

## 12. Structural Biology
The structural module is a Bio.PDB-style API for coordinate data, geometry, contact analysis, mmCIF/PDB I/O, and visualization.

* Core structure types: `Atom`, `Residue`, `Chain`, `Model`, `Structure`, `SuperpositionResult`, `AtomKDTree`
* File I/O: `read_pdb`, `read_mmcif`, `write_pdb`, `write_mmcif`
* Structural hierarchy helpers: `structure_models`, `structure_chains`, `structure_residues`, `structure_atoms`
* Coordinate access: `atom_coordinates`, `coordinate_matrix`
* Geometry and fitting: `kabsch`, `rmsd`, `superpose`, `superpose!`
* Backbone geometry: `torsion_angle`, `phi_psi`, `backbone_torsions`
* Structural metrics: `atomic_mass`, `center_of_mass`, `bounding_box`, `radius_of_gyration`, `atom_distance_matrix`, `residue_distance_matrix`, `chain_contact_matrix`, `structure_summary`, `chain_summary`
* Residue and atom selection: `select_atoms`, `select_residues`, `collapse_altlocs`, `AtomSelectionPolicy`, `residue_property`, `residue_bfactors`, `flexible_residues`, `sequence_from_structure`
* Spatial and contact analysis: `residue_contacts`, `contact_map`, `interface_residues`, `build_atom_kdtree`, `atoms_within_radius`
* Solvent accessibility and interfaces: `atom_sasa`, `residue_sasa`, `chain_sasa`, `structure_sasa`, `sasa_profile`, `residue_free_sasa`, `buried_surface_area`, `interface_profile`
* Disulfides and hydrogen bonds: `disulfide_bonds`, `HydrogenBond`, `hydrogen_bonds`
* Secondary-structure annotation: `DSSPEntry`, `read_dssp`, `annotate_dssp!`
* Backbone quality: `ramachandran_region`, `ramachandran_profile`
* Rotamers: `chi_angles`, `rotamer_label`, `rotamer_state`, `rotamer_statistics`
* Ensemble and trajectory analysis: `superpose_models!`, `ensemble_rmsd_matrix`, `trajectory_statistics`
* Contact visualization: `structure_contact_graph`, `plot_contact_graph`, `plot_contact_map`, `contact_map_svg`, `write_contact_map_svg`, `structure_contact_mermaid`, `write_structure_mermaid`
* Structure viewers and plotting: `plot_structure_atoms!`, `plot_structure_atoms`, `plot_backbone_trace!`, `plot_backbone_trace`, `plot_chain_ribbon!`, `plot_chain_ribbon`, `residue_pick_hooks`, `connect_residue_picking!`, `plot_structure_viewer`
* External tools: `run_dssp`, `run_pdb2pqr`

### Structural fidelity and performance notes
The structural parser and selection layer includes quoted-field mmCIF tokenization, preserved raw category blocks, category-specific metadata round-tripping, and parenthesized boolean selector parsing. The implementation also uses cached SASA sphere samples, reused coordinate matrices for ensemble work, and bucketed residue-contact candidate searches to reduce repeated work on large structures.

## 13. Phylogenetics and Tree Analysis
Tree-building and tree I/O cover common phylogenetic workflows.

* Tree type: `PhyloTree`
* Distances and construction: `distance_matrix`, `neighbor_joining`, `neighbor_joining_tree`, `upgma`
* Tree parsing and writing: `parse_newick`, `write_newick`, `parse_tree`, `write_tree`, `parse_phyloxml`, `write_phyloxml`, `parse_nexus`, `write_nexus`, `parse_nexml`, `write_nexml`
* Likelihood models: `JC69`, `K80`, `HKY85`, `transition_probability`, `felsenstein_likelihood`, `maximum_likelihood_tree`
* Tree traversal and analysis: `get_terminals`, `get_nonterminals`, `tree_distance`, `prune`, `root_with_outgroup`, `reroot`, `midpoint_root`, `get_parent`, `lowest_common_ancestor`, `common_ancestor`, `is_monophyletic`, `robinson_foulds_distance`, `set_metadata!`, `annotate_tree!`, `is_terminal`, `is_parent_of`, `get_path`, `trace`, `is_bifurcating`, `is_preterminal`, `total_branch_length`, `depths`, `find_clades`, `ladderize`, `count_terminals`, `collapse_clades`
* Consensus and bootstrap support: `parsimony_score`, `maximum_parsimony_tree`, `parsimony_tree`, `bootstrap_trees`, `tree_consensus`, `consensus_tree`, `strict_consensus_tree`, `majority_consensus_tree`, `bootstrap_consensus_tree`, `bootstrap_support`
* Tree rendering: `draw_ascii`, `draw_unicode`, `tree_to_dot`, `tree_to_mermaid`

## 14. Population Genetics
Population genetics includes allele statistics, LD, selection scans, kinship, simulations, and GenePop helpers.

* Core containers: `Locus`, `PopGenIndividual`, `Population`, `GenePopRecord`
* Allele and genotype summaries: `allele_frequencies`, `genotype_frequencies`, `heterozygosity_observed`, `heterozygosity_expected`
* Hardy-Weinberg and $F$-statistics: `hardy_weinberg_test`, `hardy_weinberg_exact`, `f_statistics`, `g_statistics`
* Variance partitioning and distances: `amova`, `genetic_distance`, `population_pca`, `population_pcoa`, `mantel_test`, `mismatch_distribution`
* Effective population size and migration: `estimate_ne_ld`, `estimate_ne_temporal`, `migration_rate`
* Linkage disequilibrium and phase inference: `linkage_disequilibrium`, `ld_mapping`, `ld_decay`, `infer_phase_em`
* Neutrality and selection: `segregating_sites`, `nucleotide_diversity`, `watterson_theta`, `tajimas_d`, `fu_li_d`, `fu_li_f`, `ehh`, `ihs`, `xp_ehh`, `ewens_watterson_test`, `sweepfinder_clr`
* Admixture and relatedness: `f3_statistic`, `f4_statistic`, `patterson_d`, `genetic_relationship_matrix`, `linear_mixed_model_scan`, `inbreeding_coefficient`, `relatedness`
* GenePop I/O and manipulation: `read_genepop_record`, `read_genepop`, `write_genepop`, `split_in_pops`, `split_in_loci`, `remove_population!`, `remove_locus_by_position!`, `remove_locus_by_name!`, `remove_population`, `remove_locus_by_position`, `remove_locus_by_name`
* External population tools: `run_genepop`, `run_fastsimcoal`, `run_fdist`
* Simulation helpers: `wright_fisher_simulation`, `wright_fisher_metapopulation`, `simulate_coalescent`, `site_frequency_spectrum`

## 15. GWAS and Polygenic Scoring
The GWAS layer builds directly on population genetics and adds end-to-end association workflows.

* Core data model: `GWAS`, `GenotypeMatrix`, `GWASResult`, `MetaAnalysisResult`
* PLINK IO: `read_plink`, `write_plink`
* Association scanning: `gwas_linear_scan`, `gwas_lmm_scan`
* Polygenic scoring: `calculate_prs`, `prs_ldpred`, `prs_cross_validation`
* Filtering and integration: `ld_clumping`, `meta_analyze`, `overlap_gwas_peaks`, `gene_based_test`

The module is intentionally orchestration-heavy: it keeps the matrix math in the association scanners, then reuses `GenomicRanges`, `Epigenetics`, `SystemsBio`, and `Enrichment` for downstream interpretation and module-level analysis.

## 16. Database, Classification, and Reference Modules
BioToolkit includes compact, Biopython-style helpers for common database and reference workflows.

* Restriction enzymes: `Restriction`, `RestrictionEnzyme`, `RestrictionSite`, `restriction_enzymes`, `restriction_enzyme`, `restriction_enzyme_names`, `restriction_catalog`, `restriction_sites`, `restriction_digest_map`, `find_restriction_sites`, `digest_sequence`
* Entrez / E-utilities: `Entrez`, `EntrezSearchResult`, `EntrezPostResult`, `entrez_search`, `entrez_search_ids`, `entrez_search_count`, `entrez_fetch`, `entrez_fetch_fasta`, `entrez_fetch_genbank`, `entrez_summary`, `entrez_post`, `entrez_post_ids`, `entrez_post_webenv`, `entrez_post_query_key`, `entrez_link`, `entrez_elink`, `entrez_link_ids`, `entrez_link_linksets`, `entrez_link_records`, `entrez_pubmed_search`, `entrez_pubmed_fetch`, `entrez_nuccore_fetch`, `entrez_nucleotide_search`, `entrez_protein_search`, `entrez_gene_search`, `entrez_taxonomy_search`, `entrez_genome_search`, `entrez_genome_fetch`, `parse_entrez_search_response`, `parse_entrez_post_response`
* MEDLINE / PubMed text and XML parsing: `Medline`, `MedlineRecord`, `parse_medline`, `parse_medline_xml`, `parse_medline_text`, `read_medline`
* KEGG flat-file parsing: `KEGG`, `KEGGRecord`, `KEGGPathwayRecord`, `KEGGEnzymeRecord`, `read_kegg_record`, `read_kegg_pathway`, `read_kegg_enzyme`, `kegg_field`, `kegg_entries`, `kegg_entry_id`
* Pathway graph helpers: `Pathway`, `PathwayNode`, `PathwayEdge`, `PathwayGraph`, `read_pathway_graph`, `pathway_nodes`, `pathway_edges`, `kegg_pathway_mermaid`, `write_kegg_pathway_mermaid`
* Structural classification parsing: `SCOP`, `SCOPRecord`, `read_scop_records`, `parse_scop_record`, `CATH`, `CATHRecord`, `read_cath_records`, `parse_cath_record`
* Classification hierarchy fields: `SCOPRecord.hierarchy`, `CATHRecord.hierarchy`
* Compass / EMBOSS command wrappers: `Compass`, `run_needle`, `run_water`, `run_transeq`, `run_revseq`, `run_compseq`, `run_seqret`

## 17. Plotting, Enrichment, and Single-Cell Analysis
The latest analytical modules extend BioToolkit toward an end-to-end transcriptomics and visualization stack.

* Plotting helpers: `BioPlotting`, `VolcanoPoint`, `VolcanoPlotResult`, `MAPlotResult`, `ClusteredHeatmapResult`, `ManhattanPoint`, `QQPoint`, `ForestPoint`, `ManhattanPlotResult`, `QQPlotResult`, `ForestPlotResult`, `volcano_data`, `volcano_plot`, `ma_data`, `ma_plot`, `clustered_heatmap`, `export_plot`, `manhattan_data`, `manhattan_plot`, `qq_data`, `qq_plot`, `gwas_forest_plot`
* Enrichment helpers: `Enrichment`, `IDMapper`, `EnrichmentTerm`, `EnrichmentDatabase`, `EnrichmentResult`, `load_annotation_database`, `save_annotation_database`, `build_annotation_database`, `map_id`, `map_ids`, `enrichment_test`, `go_enrichment`, `kegg_enrichment`, `dotplot`
* Enrichment fixture helpers: `builtin_annotation_terms`, `builtin_annotation_database`
* Single-cell helpers: `SingleCell`, `SingleCellExperiment`, `count_matrix`, `normalize_counts`, `sctransform`, `run_pca`, `run_umap`, `cluster_cells`, `find_cluster_markers`, `summarize_clusters`, `cluster_marker_summary`

## 18. Epigenomics Engine
The `Epigenetics` module adds a unified object model and a native epigenomics workflow layer on top of the existing BioToolkit interval, RNA-seq, single-cell, motif, and HMM primitives.

* Core data model: `Epigenetics`, `SparseCoverageVector`, `Peak`, `PeakSet`, `Epigenome`
* Coverage and peaks: `calculate_coverage`, `coverage_depth`, `coverage_segments`, `call_peaks`, `normalize_gc_bias`
* Differential regulation: `count_overlaps`, `differential_binding`
* Methylation: `MethylationCall`, `MethylationExperiment`, `MethylationResult`, `bin_methylation`, `differential_methylation`
* Single-cell chromatin: `SingleCellChromatinExperiment`, `tfidf`, `run_lsi`, `rsvd`, `gene_activity_score`, `calculate_coaccessibility`
* Motifs and footprints: `compute_motif_deviations`, `detect_footprints`
* 3D genomics: `ContactMatrix`, `directionality_index`, `detect_tads`

The module is designed to be layered rather than standalone:

* `Epigenome` stores intervals, coverage, counts, and sample metadata together.
* Coverage is computed with a sweep-line representation instead of per-base iteration.
* Peak calling, differential binding, and methylation reuse the existing sparse statistical stack.
* Single-cell chromatin and Hi-C helpers reuse the package’s sparse matrix and HMM infrastructure.

## 19. Overall Design
BioToolkit keeps the public API broad but consistent:

* Core data structures are exported directly for interactive use.
* Parsing and writing functions are provided for the most common bioinformatics file formats.
* Structural biology features cover PDB and mmCIF round-tripping, contact analysis, solvent accessibility, and visualization helpers.
* Phylogenetics and population genetics are available in the same namespace for end-to-end workflow scripting.
* Database/reference workflows now include restriction enzymes, Entrez, PubMed/MEDLINE, KEGG, SCOP, CATH, and Compass-style wrappers.

For implementation details, see [src/BioToolkitBody.jl](src/BioToolkitBody.jl) and the repository-level package note in [src/BioToolkit.md](src/BioToolkit.md).

## 20. Benchmark and Parity Coverage
BioToolkit is accompanied by a comparison suite that keeps the public API honest against Biopython where a direct baseline exists.

* Sequence toolkit coverage: [scripts/sequence_benchmark.jl](scripts/sequence_benchmark.jl), [Examples/Comparison/sequence_toolkit_compare.jl](Examples/Comparison/sequence_toolkit_compare.jl), [Examples/Comparison/hamming_compare.jl](Examples/Comparison/hamming_compare.jl), [Examples/Comparison/expanded_biopython_compare.jl](Examples/Comparison/expanded_biopython_compare.jl). The measured overlap now includes DNA transcription, melting temperature, DNA molecular weight, and sequence slicing in addition to count/GC/reverse-complement/translation/hamming/k-mer timing.
* Alignment and scoring coverage: [Examples/Comparison/pairwise_align_compare.jl](Examples/Comparison/pairwise_align_compare.jl), [Examples/Comparison/pairwise_affine_compare.jl](Examples/Comparison/pairwise_affine_compare.jl), [Examples/Comparison/protein_matrix_compare.jl](Examples/Comparison/protein_matrix_compare.jl), [Examples/Comparison/local_align_compare.jl](Examples/Comparison/local_align_compare.jl)
* File-format coverage: [Examples/Comparison/fastq_compare.jl](Examples/Comparison/fastq_compare.jl), [Examples/Comparison/fasta_index_compare.jl](Examples/Comparison/fasta_index_compare.jl), [Examples/Comparison/genbank_compare.jl](Examples/Comparison/genbank_compare.jl), [Examples/Comparison/genbank_write_compare.jl](Examples/Comparison/genbank_write_compare.jl), [Examples/Comparison/annotation_compare.jl](Examples/Comparison/annotation_compare.jl)
* Quality coverage: [Examples/Comparison/quality_compare.jl](Examples/Comparison/quality_compare.jl) includes Phred decode, mean quality, trimming, and Phred string encoding.
* Flat-file parsing coverage: [Examples/Comparison/flatfile_compare.jl](Examples/Comparison/flatfile_compare.jl) adds a direct EMBL parse comparison against Biopython on a synthetic record.
* Trace parsing coverage: [Examples/Comparison/abif_compare.jl](Examples/Comparison/abif_compare.jl) times the optimized ABIF parser on a synthetic trace fixture.
* Search and parsing coverage: [Examples/Comparison/expanded_biopython_compare.jl](Examples/Comparison/expanded_biopython_compare.jl) adds BLAST XML/tabular parsing, sequence slicing, and sequence utilities that were missing from the earlier comparison set.
* Motif and MSA coverage: [Examples/Comparison/motif_pwm_compare.jl](Examples/Comparison/motif_pwm_compare.jl), [Examples/Comparison/motif_scan_compare.jl](Examples/Comparison/motif_scan_compare.jl), [Examples/Comparison/msa_compare.jl](Examples/Comparison/msa_compare.jl), [Examples/Comparison/msf_compare.jl](Examples/Comparison/msf_compare.jl), [Examples/Comparison/clustal_muscle_compare.jl](Examples/Comparison/clustal_muscle_compare.jl)
* Phylogenetics coverage: [Examples/Comparison/nj_compare.jl](Examples/Comparison/nj_compare.jl), [Examples/Comparison/midpoint_root_compare.jl](Examples/Comparison/midpoint_root_compare.jl), [Examples/Comparison/consensus_support_newick_compare.jl](Examples/Comparison/consensus_support_newick_compare.jl), [Examples/Comparison/prune_reroot_compare.jl](Examples/Comparison/prune_reroot_compare.jl), [Examples/Comparison/graph_export_compare.jl](Examples/Comparison/graph_export_compare.jl), [Examples/Comparison/phylip_compare.jl](Examples/Comparison/phylip_compare.jl), [Examples/Comparison/expanded_julia_only_compare.jl](Examples/Comparison/expanded_julia_only_compare.jl)
* Structural coverage: [Examples/Comparison/expanded_biopython_compare.jl](Examples/Comparison/expanded_biopython_compare.jl) includes the 10,000-atom neighbor-search and superposition checks, while the Julia-only companion keeps larger structural and phylogenetic runs measurable even when the Python baseline is optional.
* Genomic coverage: [Examples/Comparison/window_coverage_compare.jl](Examples/Comparison/window_coverage_compare.jl), [Examples/Comparison/cuda_histogram_compare.jl](Examples/Comparison/cuda_histogram_compare.jl). CUDA-backed helpers are now intended to load on demand instead of during the base package startup path.
* Database and reference workflows: [benchmark_report.md](benchmark_report.md) summarizes the restriction, Entrez, MEDLINE, KEGG, SCOP, CATH, GenBank, and annotation overlap checks and separates them from parser-first Julia-only surfaces.
* Protein coverage: [Examples/Comparison/protein_stats_compare.jl](Examples/Comparison/protein_stats_compare.jl) covers the individual ProtParam-style calculations plus the bundled `protparam` workflow.
* Population genetics coverage: [Examples/Comparison/popgen_compare.jl](Examples/Comparison/popgen_compare.jl) and [scripts/bench_popgen.jl](scripts/bench_popgen.jl)

The goal of the comparison suite is not to force every feature into a Biopython mold. It is to keep the shared surface measured where a meaningful Python baseline exists, while still documenting the Julia-only analytical paths separately in [benchmark_report.md](benchmark_report.md).

The current design deliberately keeps the root package body thin and pushes feature-specific logic into the nested `BioToolkit.jl/src/` modules. Optional CUDA entrypoints now hang off a dedicated lazy loader, so the interval, RNA-seq, structure, motif, population-genetics, and GPU paths can evolve independently without adding startup-time imports.