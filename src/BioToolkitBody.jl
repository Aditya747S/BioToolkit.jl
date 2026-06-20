using Arrow
using Tables
using LinearAlgebra
using Distributions
using Statistics
using Random
using SpecialFunctions
using Printf
using BGZFStreams
using BigWig

include("analysisresults.jl")
include("biotypes.jl")
include("schema.jl")
include("io.jl")
include("lazy_gpu.jl")
include("query.jl")
include("genomicranges.jl")
include("bam.jl")
include("record.jl")
include("quality.jl")
include("sequence.jl")
include("longread.jl")
include("protein.jl")
include("annotation.jl")
include("align.jl")
include("msa.jl")
include("motif.jl")
include("search.jl")
include("databases.jl")
include("hmm.jl")
include("phylo.jl")
include("structure.jl")
include("coevolution.jl")
include("gene_prediction.jl")
include("popgen.jl")
include("differentialexpression.jl")
include("enrichment.jl")
include("singlecell.jl")
include("immunology.jl")
include("spatial.jl")
include("deeplearning.jl")
include("epigenetics.jl")
include("clinical.jl")
include("microbiome.jl")
include("proteomics.jl")
include("metabolomics.jl")
include("systemsbio.jl")
include("somatic.jl")
include("crispr.jl")
include("gwas.jl")
include("bioplotting.jl")
include("bioconductor_compat.jl")
include("pipeline.jl")

using .Restriction: RestrictionEnzyme, RestrictionSite, restriction_enzymes, restriction_enzyme, restriction_enzyme_names, restriction_catalog, restriction_sites, restriction_digest_map, find_restriction_sites, digest_sequence
using .Entrez: EntrezSearchResult, EntrezPostResult, entrez_search, entrez_search_ids, entrez_search_count, entrez_fetch, entrez_fetch_fasta, entrez_fetch_genbank, entrez_summary, entrez_post, entrez_post_ids, entrez_post_webenv, entrez_post_query_key, entrez_link, entrez_elink, entrez_link_ids, entrez_link_linksets, entrez_link_records, entrez_pubmed_search, entrez_pubmed_fetch, entrez_nuccore_fetch, entrez_nucleotide_search, entrez_protein_search, entrez_gene_search, entrez_taxonomy_search, entrez_genome_search, entrez_genome_fetch, parse_entrez_search_response, parse_entrez_post_response
using .Medline: MedlineRecord, parse_medline, parse_medline_xml, parse_medline_text, read_medline
using .KEGG: KEGGRecord, KEGGPathwayRecord, KEGGEnzymeRecord, read_kegg_record, read_kegg_pathway, read_kegg_enzyme, kegg_field, kegg_entries, kegg_entry_id
using .Pathway: PathwayNode, PathwayEdge, PathwayGraph, read_pathway_graph, pathway_nodes, pathway_edges, kegg_pathway_mermaid, write_kegg_pathway_mermaid
using .GenomicRanges: GenomicRanges, GenomicInterval, IntervalCollection, CoverageSegment, build_collection, read_intervals, overlap, find_overlaps, find_overlaps_parallel, nearest, find_nearest, follow, precede, coverage, complement
using .LongRead: NanoporeRead, PacBioRead, StructuralVariant, MinimizerSeed, CandidateRegion, BandedAlignmentResult, canonical_kmer_hash, minimizer_sketch, build_minimizer_index, find_candidate_regions, banded_semiglobal_alignment, split_read_evidence, discordant_pair_evidence, cluster_structural_variants, call_structural_variants
using .DifferentialExpression: DifferentialExpression, CountMatrix, DEResult, GLMSolver, DESeqDataSet, DESeqDataSetFromMatrix, counts, design, calc_tmm_factors, calc_norm_factors, estimateSizeFactorsForMatrix, estimate_dispersions, estimate_dispersions_prior, cooks_distance, replace_outliers, estimateSizeFactors, estimateDispersions, estimateDispersionsGeneEst, estimateDispersionsFit, estimateDispersionsMAP, estimateDispersionsPriorVar, nbinomWaldTest, nbinomLRT, estimate_dispersion_edgeR, trend_dispersion_edgeR, dispersion_prior_dof, exact_test_edgeR, glm_ql_fit, glm_ql_f_test, edgeR_qlf_test, DESeq, results, resultsNames, sizeFactors, sizeFactors!, normalizationFactors, normalizationFactors!, makeExampleDESeqDataSet, benjamini_hochberg, differential_expression, filter_low_counts, shrink_lfc, vst, fit_gene_fast!, remove_batch_effect, combat_correction, estimate_surrogates
using .Enrichment: IDMapper, EnrichmentTerm, EnrichmentDatabase, EnrichmentResult, load_annotation_database, save_annotation_database, build_annotation_database, builtin_annotation_database, builtin_annotation_terms, map_id, map_ids, enrichment_test, go_enrichment, kegg_enrichment, dotplot
using .SingleCell: SingleCellExperiment, SingleCellProjectionModel, count_matrix, normalize_counts, sctransform, find_variable_features, fit_singlecell_projection_model, project_singlecell, run_pca, find_neighbors, find_spatial_neighbors, run_umap, cluster_cells, find_clusters, calculate_pseudotime, find_cluster_markers, find_markers, summarize_clusters, cluster_marker_summary, detect_doublets, integrate_batches, integrate_data, score_cell_cycle, SpatialMoranResult, attach_spatial_coords!, moran_i_test, find_spatially_variable_genes, RNAVelocityResult, DynamicalRNAVelocityResult, attach_velocity_layers!, calculate_rna_velocity, calculate_dynamical_rna_velocity, plot_trajectory, plot_rna_velocity, plot_rna_velocity_quiver, plot_communication_network, plot_communication_pathways, plot_ligand_receptor_report, gpu_singlecell_experiment, materialize_singlecell_experiment, save_singlecell_experiment, load_singlecell_experiment, SingleCellArchive, save_singlecell_archive, load_singlecell_archive, archive_to_singlecell_experiment, SingleCellArchiveBrowser, browse_singlecell_archive, archive_layer, archive_embedding, LigandReceptorPair, CellCommunicationResult, CellCommunicationNetwork, CommunicationPathwaySummary, LigandReceptorReport, default_ligand_receptor_pairs, find_cell_communication, communication_network, communication_pathway_summary, rank_ligand_receptor_report, WNNResult, weighted_nearest_neighbors, PerturbationPredictionResult, predict_perturbation, read_h5ad, write_h5ad, CellTypeAnnotationResult, annotate_cell_types, AmbientRNARemovalResult, remove_ambient_rna, SingleCellViewer, interactive_singlecell_viewer, lasso_select_cells, recluster_singlecell_viewer!, cell_hover_text, top_expressed_genes
using .SpatialDeconvolution: SpatialExperiment, DeconvolutionResult, build_reference_matrix, rctd_deconvolution, cell2location_deconvolution, spatially_variable_genes, spatial_autocorrelation, spatial_neighborhood_graph, spatial_domain_clustering, ligand_receptor_spatial_score, build_spot_deconvolution_qc, spatial_coexpression_modules, mark_tissue_boundary_spots, spatial_pseudotime
using .Epigenetics: SparseCoverageVector, Peak, PeakSet, PeakSupport, Epigenome, MethylationCall, MethylationExperiment, MethylationResult, SingleCellChromatinExperiment, CoaccessibilityEdge, TadResult, ContactMatrix, calculate_coverage, coverage_depth, coverage_segments, call_peaks, normalize_gc_bias, summarize_peak_support, count_overlaps, differential_binding, bin_methylation, differential_methylation, tfidf, run_lsi, rsvd, gene_activity_score, calculate_coaccessibility, compute_motif_deviations, detect_footprints, directionality_index, detect_tads
using .Clinical: PatientCohort, KaplanMeierResult, CoxResult, CoxTermResult, MAFRecord, MAFSummary, ROCResult, CIFResult, NeuralCoxResult, DoseResponseResult, OncoprintResult, read_maf, summarize_maf, tcga_query, tcga_download_files, merge_tcga_count_files, tcga_ingest, kaplan_meier, kaplan_meier_plot, logrank_test, cox_ph, forest_plot, survival_roc, cif_curve, neural_cox, dose_response_curve, oncoprint
using .Microbiome: CommunityProfile, PCoAResult, NMDSResult, ANCOMResult, SongbirdResult, MicrobiomeNetwork, SourceTrackingResult, clr_transform, ilr_transform, bray_curtis, unifrac, weighted_unifrac, pairwise_bray_curtis, pairwise_unifrac, shannon_entropy, simpson_index, faith_pd, pcoa, pcoa_plot, nmds, ancom, songbird, cooccurrence_network, network_plot, source_tracking_model, source_tracking, source_tracking_posterior_summary
using .Proteomics: MassSpecExperiment, Spectrum, MassSpecPeak, PeakDetectionResult, AlignmentResult, DifferentialAbundanceResult, SparsePLSDAResult, DeNovoResult, read_mzml, detect_peaks, align_samples, qrilc_impute, differential_abundance, mixed_model_abundance, sparse_pls_da, stream_mass_spec, build_spectrum_graph, de_novo_sequence
using .Metabolomics: MetabolomicsSourceTrackingResult, metabolomics_source_tracking_model, metabolomics_source_tracking, metabolomics_source_tracking_posterior_summary, metabolomics_differential_abundance, metabolomics_streaming_analysis, annotate_metabolite_features, quantify_metabolite_variation, MetaboliteAnnotation, IsotopeTrace, MetabolicPathwayResult, normalize_metabolomics, annotate_adducts, isotope_tracer_analysis, metabolite_pathway_enrichment, quality_control_metabolomics, pca_metabolomics, metabolite_correlation_network, fold_change_metabolomics, volcano_metabolomics, missing_value_analysis, metabolite_classification, nmr_peak_table, batch_effect_correction_metabolomics, targeted_quantification, metabolomics_power_analysis, feature_clustering_metabolomics
using .SystemsBio: GeneNetwork, SoftThresholdResult, LimmaResult, VoomResult, NetworkInferenceResult, MultiOmicsFactorAnalysisResult, pick_soft_threshold, find_modules, module_eigengenes, module_dendrogram, voom_transform, voom_with_quality_weights, estimate_limma_hyperparameters, limma_moderate_ttest, limma_fit, limma_deresults, eBayes, contrasts_fit, duplicateCorrelation, remove_batch_effect_limma, gsea, infer_network, multi_omics_factor_analysis
using .GWAS: GenotypeMatrix, GWASResult, MetaAnalysisResult, BedReader, BgenReader, BgenVariant, read_plink, write_plink, read_bed_genotypes, read_bgen, write_bgen, gwas_linear_scan, gwas_lmm_scan, gwas_logistic_scan, gwas_gxe_interaction, calculate_maf, calculate_missingness, missingness_filter, info_score_filter, gwas_qc_report, hwe_exact, calculate_hwe_pvalues, hwe_filter, genomic_control_lambda, apply_genomic_control, calculate_ld_matrix, calculate_grm, loco_projection, calculate_prs, prs_ldpred, prs_cross_validation, ld_clumping, meta_analyze, overlap_gwas_peaks, gene_based_test, to_plink_dataframe, write_plink_sumstats, gwas_migration_guide, LDSCResult, ldsc_heritability, estimate_heritability_greml, partitioned_heritability, ldsc_genetic_correlation, SuSiEResult, fine_map_susie, calculate_credible_set, posterior_inclusion_probability, MRResult, mr_two_sample, mr_egger, mr_pleiotropy_test, conditional_analysis, cojo_stepwise, joint_analysis, gwas_pca, project_pca, calculate_ibd, calculate_king_kinship, detect_related_pairs, coloc_result, coloc_abf, ihs_score, xp_ehh_score, fst_outlier_test, rank_inverse_normal, rank_inverse_normal!, calculate_sex_check, calculate_heterozygosity_outliers, sample_qc_report, normalise_chromosome, filter_variants, filter_samples, merge_genotype_matrices, read_oxford_gen, compute_ld_scores, prune_ld, flip_alleles, harmonise_alleles, sample_prune_relatedness, ancestry_inference, locus_zoom, score_test_linear, likelihood_ratio_test, gwas_survival_scan, gwas_ordinal_scan, gwas_multivariate_scan, burden_test, skat_test, conditional_fdr, storey_pi0_estimate, gwas_power_calculation, permutation_pvalue, genomic_sem_fit, twas_scan, smr_test, liftover, functional_annotation, ld_expand_credible_set, estimate_effective_n, he_regression_variance_components, reml_variance_components, simulate_genotypes, simulate_phenotype, dosage_to_hardcall, info_score_from_dosage, popcorn_genetic_correlation, ebi_lookup, save_gwas_result, load_gwas_result
using .BioconductorCompat: biocparallel_map, scater_qc_metrics, scran_pool_size_factors, scuttle_aggregate_counts, singler_annotate, scmap_project, splatter_simulate_counts, diffbind_differential_binding, csaw_window_differential_binding, chromvar_deviations, minfi_dmp, bsseq_dmr, xcms_peak_workflow, lipidr_differential_abundance, mixomics_factor_analysis, mofa2_factor_analysis, variantannotation_annotate, varianttools_filter_variants, gwascat_lookup, complexheatmap_payload, gviz_track_table
using .SingleCell: bayesspace_like_domains, spagc_like_domains, sparkx_spatial_de, spatialde_gp_de, spatial_trajectory_graph, niche_weighted_communication, cell2location_like_segmentation, mixscape_like_contrast, spatial_markov_refine_domains, spatial_lr_permutation_test, mixscape_multibatch_contrast, perturbseq_pseudobulk, milo_like_neighborhood_da, perturbation_synergy_scores
using .Immunology: extract_cdr3, clonotype_table, assign_vdj_segments, isotype_switching_summary, germline_usage_bias, bepipred_like_scores, mhcflurry_like_scores, repertoire_diversity_metrics, clonal_expansion_test, cdr3_length_spectrum, tcrdist_like_matrix, igblast_assign, imgt_highvquest_manifest, paired_clonotype_table, differential_vj_usage, somatic_hypermutation_rate, bcr_affinity_maturation_score, antigen_specificity_clustering, clonotype_trajectory, v_gene_family_usage, hla_type_enrichment, cdr3_physicochemical_properties, network_centrality_clonotypes, convergent_clonotype_detection, lineage_tree_from_clones, epitope_binding_profile
using .SystemsBio: mofa_plus_integration, cca_integration, anchor_based_integration, totalvi_like_integration, causal_grn_inference, gene_program_nmf, mofa_plus_em, sparse_cca_integration, causal_ate_regression
using .LongRead: isoseq_consensus, minimap2_align, graphaligner_align, sniffles_call, svim_call, whatshap_phase, phase_reads_by_alleles, overlap_layout_consensus, read_n50, overlap_graph_table, porec_contact_table, omnic_contact_matrix
using .DeepLearning: scvi_like_embedding, cellassign_like_mapping, geneformer_like_embedding, scgpt_like_embedding, scbert_like_embedding, attention_grn, batch_corrected_latent, contrastive_cell_embedding, flux_autoencoder_embedding, flux_mlp_classifier, scgen_like_perturbation, graphsca_label_transfer, cell_type_denoising, sparse_autoencoder_features, trajectory_neural_ode, multimodal_wnn_embedding, protein_sequence_embedding, zero_shot_cell_annotation, gene_regulatory_network_gnn, self_supervised_pretraining, cell_cycle_regression, deep_factorization_embedding
using .Epigenetics: chromhmm_like_segmentation, hic_ice_normalize, hic_kr_normalize, hic_ab_compartments, tobias_like_footprints, parse_bismark_coverage, insulation_score, compartment_boundary_candidates
using .Proteomics: dia_like_quantification, phosphosite_localization, glycoproteomics_motif_table, project_peptides_to_structure, protein_inference_top3, ptm_site_enrichment
using .Clinical: pharmacogenomics_star_alleles, omop_visit_summary, synthpop_like_cohort, trial_suitability_scores, cpic_metabolizer_phenotype, pharmgkb_like_recommendations, propensity_score_match
using .Microbiome: mag_bin_contigs, kraken_like_classify, viral_contig_scores, humann_like_pathways, strainge_like_variants, lca_taxonomy_from_votes, strain_haplotype_profile
using .Pipeline: terra_workspace_plan, anvil_workspace_manifest, distributed_pipeline_map, retry_execute_pipeline, containerized_node, slurm_array_plan
using .Somatic: germline_bayesian_call, somatic_tumor_normal_call, haplotypecaller_call, mutect2_call, strelka2_call, sv_breakpoint_graph_table, vep_like_annotation, tumor_mutational_burden, mutational_signature_nmf, cbs_like_segments, cosmic_signature_attribution, clonal_deconvolution_pyclone, clonal_evolution_tree, allelic_imbalance_test, copy_number_from_coverage, sv_fusion_candidates, driver_enrichment, variant_tier_classification, somatic_hotspot_scan, mutational_spectrum_96
using .BioPlotting: VolcanoPoint, VolcanoPlotResult, MAPoint, MAPlotResult, ClusteredHeatmapResult, ManhattanPoint, ManhattanPlotResult, QQPoint, QQPlotResult, ForestPoint, ForestPlotResult, volcano_data, volcano_plot, ma_data, ma_plot, clustered_heatmap, export_plot, manhattan_data, manhattan_plot, qq_data, qq_plot, gwas_forest_plot
using .Pipeline: FileArtifact, PipelineNode, PipelineGraph, file_artifact, pipeline_node, add_node!, add_dependency!, validate_pipeline, execution_levels, execute_pipeline, align_reads, count_features, template_read_fasta_node, template_align_reads_node, template_count_features_node, template_differential_expression_node
using .Enrichment: IDMapper, EnrichmentTerm, EnrichmentDatabase, EnrichmentResult, load_annotation_database, save_annotation_database, build_annotation_database, builtin_annotation_database, builtin_annotation_terms, map_id, map_ids, enrichment_test, go_enrichment, kegg_enrichment, dotplot, GSEAResult, fgsea_like, gsea_preranked, GSVAResult, gsva_score, network_propagation, heat_diffusion_enrichment, reactome_enrichment, wikipathways_enrichment, msigdb_enrichment, gene_set_overlap_matrix, jaccard_similarity_matrix, leading_edge_genes, enrichment_map, competitive_gene_set_test, self_contained_gene_set_test, enrichment_heatmap_data, bubble_chart_data, gene_ontology_semantic_similarity, go_slim_mapping, term_to_gene_matrix, rank_genes_by_set_membership
using .Coevolution: ContactMap, PseudoLikelihoodModel, filter_alignment_for_dca, sequence_reweighting, fit_pseudolikelihood_model, compute_contact_scores, predict_contact_map, top_contact_pairs, fold_from_contacts, mutual_information_contacts, direct_information_contacts, column_conservation_scores, sequence_logo_entropy, evolutionary_coupling_network, contact_enrichment_statistics, shrinkage_precision_contacts, phylogenetic_correction, positional_covariation_matrix, gap_analysis, contact_precision_recall, alignment_quality_report
using .SCOP: SCOPRecord, read_scop_records, parse_scop_record
using .CATH: CATHRecord, read_cath_records, parse_cath_record
using .CRISPR: GuideRNA, CRISPRSystem, OffTarget, EditingWindow, SpCas9, SaCas9, Cas12a, CasX, SpRY, BE3, ABE8e, CBE4max, design_guides, score_on_target, find_pam_sites, enumerate_off_targets, cfd_score, mit_score, guide_gc_content, design_base_editor_guides, analyze_editing_window, design_pegrna, prime_editing_guide_score, crispr_screen_analysis, mageck_like_test, design_library, library_coverage_stats, design_hdr_template, predict_indels, indel_distribution, filter_guides, rank_guides, guide_specificity_score, genome_wide_offtarget_summary, pam_matches

include("browser.jl")

function _export_public_bindings!(mod::Module)
	for name in names(mod; all=true, imported=true)
		startswith(string(name), "_") && continue
		isdefined(mod, name) || continue
		value = getfield(mod, name)
		value isa Union{Function,Type} || continue
		parent = try
			parentmodule(value)
		catch
			nothing
		end
		parent === mod || continue
		@eval export $(name)
	end
end

_export_public_bindings!(@__MODULE__)
for name in names(@__MODULE__; all=true, imported=true)
	isdefined(@__MODULE__, name) || continue
	value = getfield(@__MODULE__, name)
	value isa Module || continue
	startswith(string(value), "BioToolkit.") || continue
	if value === @__MODULE__
		continue
	end
	_export_public_bindings!(value)
end
