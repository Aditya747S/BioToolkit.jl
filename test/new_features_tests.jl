using Test
using Random
using Statistics
using LinearAlgebra
using BioToolkit
using DataFrames

# =============================================================================
# Helpers shared across test groups
# =============================================================================

function _make_sce(n_genes=6, n_cells=8; seed=1)
    rng = MersenneTwister(seed)
    counts = abs.(rand(rng, Int, n_genes, n_cells) .% 100) .+ 1
    genes  = ["G$i" for i in 1:n_genes]
    cells  = ["C$i" for i in 1:n_cells]
    SingleCellExperiment(counts, genes, cells)
end

function _make_spatial(n_genes=6, n_spots=9; seed=2)
    rng = MersenneTwister(seed)
    counts = abs.(rand(rng, Int, n_genes, n_spots) .% 80) .+ 5
    genes  = ["G$i" for i in 1:n_genes]
    spots  = ["S$i" for i in 1:n_spots]
    coords = hcat(rand(rng, n_spots), rand(rng, n_spots))
    sce    = SingleCellExperiment(counts, genes, spots; spatial_coords=coords)
    SpatialExperiment(sce)
end

function _make_msa(n_seqs=40, l=12; seed=3)
    rng = MersenneTwister(seed)
    aa  = collect("ACDEFGHIKLMNPQRSTVWY")
    seqs = String[]
    for _ in 1:n_seqs
        state = rand(rng, Bool)
        chars = [rand(rng, aa) for _ in 1:l]
        chars[1] = state ? 'A' : 'G'
        chars[7] = state ? 'V' : 'L'
        push!(seqs, String(chars))
    end
    MultipleSequenceAlignment(seqs)
end

# =============================================================================
# SPATIAL TRANSCRIPTOMICS  (spatial.jl new features)
# =============================================================================
@testset "Spatial — Spatially Variable Genes (Moran's I)" begin
    sp = _make_spatial(8, 12)
    svg = spatially_variable_genes(sp; k=3, permutations=9, min_expr_frac=0.0)
    @test svg isa DataFrame
    @test hasproperty(svg, :gene_id)
    @test hasproperty(svg, :morans_i)
    @test hasproperty(svg, :z_score)
    @test hasproperty(svg, :pvalue)
    @test nrow(svg) > 0
    # Moran's I bounded between -1 and 1 (approximately)
    @test all(isfinite, svg.morans_i)
    @test all(pv -> 0.0 <= pv <= 1.0, svg.pvalue)
    # top_n filtering
    svg_top5 = spatially_variable_genes(sp; k=3, permutations=0, min_expr_frac=0.0, top_n=5)
    @test nrow(svg_top5) <= 5
end

@testset "Spatial — Spatial Autocorrelation (Moran + Geary)" begin
    sp = _make_spatial(6, 10)
    ac = spatial_autocorrelation(sp; k=3, min_expr_frac=0.0)
    @test ac isa DataFrame
    @test hasproperty(ac, :morans_i)
    @test hasproperty(ac, :gearys_c)
    @test all(isfinite, ac.morans_i)
    @test all(isfinite, ac.gearys_c)
end

@testset "Spatial — Neighbourhood Graph" begin
    sp = _make_spatial(4, 6)
    res = spatial_neighborhood_graph(sp; k=2)
    @test res isa NamedTuple
    @test res.edges isa DataFrame
    @test res.adjacency isa Matrix
    @test nrow(res.edges) > 0
    @test size(res.adjacency, 1) == size(res.adjacency, 2) == 6
    # Every spot should have at least 1 neighbour
    @test all(sum(res.adjacency, dims=2) .> 0)
end

@testset "Spatial — Domain Clustering" begin
    sp = _make_spatial(6, 12)
    dom = spatial_domain_clustering(sp; n_domains=3, k=3, seed=42)
    @test dom isa DataFrame
    @test hasproperty(dom, :spot_id)
    @test hasproperty(dom, :domain)
    @test nrow(dom) == 12
    # Labels should be integers in [1, n_domains]
    @test all(x -> 1 <= x <= 3, dom.domain)
end

@testset "Spatial — Ligand–Receptor Scoring" begin
    sp = _make_spatial(6, 8)
    lr_pairs = [("G1", "G2"), ("G3", "G4")]
    scores = ligand_receptor_spatial_score(sp, lr_pairs; k=3)
    @test scores isa DataFrame
    @test hasproperty(scores, :lr_pair)
    @test hasproperty(scores, :lr_spatial_score)
    @test nrow(scores) > 0
    @test all(isfinite, scores.lr_spatial_score)
end

@testset "Spatial — Deconvolution QC" begin
    rng = MersenneTwister(7)
    n_spots, n_types = 10, 3
    fracs = rand(rng, n_spots, n_types)
    fracs ./= sum(fracs, dims=2)
    residuals = rand(rng, n_spots) .* 5
    result = DeconvolutionResult(
        ["S$i" for i in 1:n_spots],
        ["T1","T2","T3"],
        fracs, residuals, :rctd
    )
    qc = build_spot_deconvolution_qc(result)
    @test qc isa DataFrame
    @test hasproperty(qc, :confidence)
    @test hasproperty(qc, :entropy)
    @test nrow(qc) == n_spots
    @test all(c -> c in ("singlet","mixed","low_confidence"), qc.confidence)
end

@testset "Spatial — Tissue Boundary Spots" begin
    sp = _make_spatial(4, 9)
    bnd = mark_tissue_boundary_spots(sp; k=3, threshold=0.5)
    @test bnd isa DataFrame
    @test hasproperty(bnd, :is_boundary)
    @test nrow(bnd) == 9
    @test bnd.is_boundary isa BitVector || eltype(bnd.is_boundary) == Bool
end

@testset "Spatial — Spatial Pseudotime" begin
    sp = _make_spatial(4, 8)
    pt = spatial_pseudotime(sp; root_spot=1, k=3, n_diffusion_steps=5)
    @test pt isa DataFrame
    @test hasproperty(pt, :pseudotime)
    @test nrow(pt) == 8
    @test all(isfinite, pt.pseudotime)
    @test pt.pseudotime[1] ≈ 0.0 atol=0.01    # root should be ~0
end

@testset "Spatial — Co-expression Modules" begin
    sp = _make_spatial(8, 12)
    res = spatial_coexpression_modules(sp; n_modules=3, k=3, n_svgs=8, min_expr_frac=0.0, seed=1)
    @test res isa NamedTuple
    @test res.modules isa DataFrame
    @test res.module_scores isa Matrix
    @test hasproperty(res.modules, :gene_id)
    @test hasproperty(res.modules, :module)
    @test size(res.module_scores, 1) == 12   # spots
end

# =============================================================================
# SOMATIC VARIANTS  (somatic.jl new features)
# =============================================================================
@testset "Somatic — COSMIC Signature Attribution" begin
    rng = MersenneTwister(11)
    n_ctx   = 96
    n_sigs  = 4
    n_samp  = 3
    sigs    = rand(rng, n_ctx, n_sigs)
    sigs ./= sum(sigs, dims=1)
    exposures_true = rand(rng, n_sigs, n_samp)
    exposures_true ./= sum(exposures_true, dims=1)
    catalog = sigs * exposures_true  .* 1000
    catalog = round.(Int, catalog)

    res = cosmic_signature_attribution(Float64.(catalog), sigs)
    @test res isa DataFrame
    @test hasproperty(res, :sample_index)
    @test hasproperty(res, :reconstruction_cosine)
    @test nrow(res) == n_samp
    @test all(x -> 0.0 <= x <= 1.0, res.reconstruction_cosine)

    # Vector single-sample overload
    res1 = cosmic_signature_attribution(Float64.(catalog[:, 1]), sigs)
    @test nrow(res1) == 1
end

@testset "Somatic — PyClone-like CCF EM" begin
    rng = MersenneTwister(13)
    vafs = vcat(rand(rng, 20) .* 0.2 .+ 0.05, rand(rng, 20) .* 0.15 .+ 0.35)
    res  = clonal_deconvolution_pyclone(vafs; n_clones=2, purity=0.8, n_iter=100, seed=1)
    @test res isa NamedTuple
    @test length(res.ccf_means) == 2
    @test length(res.mixing_weights) == 2
    @test sum(res.mixing_weights) ≈ 1.0 atol=1e-4
    @test isfinite(res.bic)
    @test size(res.responsibilities) == (40, 2)

    cn_vafs = [0.25, 1 / 6]
    cn_copy_numbers = [2.0, 4.0]
    cn_res = clonal_deconvolution_pyclone(cn_vafs; n_clones=1, purity=0.5, copy_numbers=cn_copy_numbers, n_iter=50, seed=2)
    @test cn_res.ccf_means[1] ≈ 1.0 atol=0.05
end

@testset "Somatic — Clonal Evolution Tree" begin
    ccf_df = DataFrame(
        sample = repeat(["S1","S2"], inner=3),
        clone  = repeat(["C1","C2","C3"], outer=2),
        ccf    = [0.9, 0.5, 0.3, 0.85, 0.45, 0.25],
    )
    res = clonal_evolution_tree(ccf_df)
    @test res isa NamedTuple
    @test res.edges isa DataFrame
    @test hasproperty(res.edges, :parent_clone)
    @test hasproperty(res.edges, :child_clone)
    @test nrow(res.edges) == 2   # 3 clones → 2 edges
end

@testset "Somatic — Allelic Imbalance Test" begin
    rng = MersenneTwister(17)
    n_loci = 50
    ref    = rand(rng, 20:60, n_loci)
    # Balanced (no LOH)
    alt_bal = rand(rng, 15:45, n_loci)

    loci_df = allelic_imbalance_test(ref, alt_bal)
    @test loci_df isa DataFrame
    @test hasproperty(loci_df, :baf)
    @test hasproperty(loci_df, :pvalue)
    @test nrow(loci_df) == n_loci
    @test all(0.0 .<= loci_df.baf .<= 1.0)

    # With segments
    seg_df = DataFrame(segment_start=[1,26], segment_end=[25,50])
    res = allelic_imbalance_test(ref, alt_bal; segments=seg_df)
    @test res isa NamedTuple
    @test res.loci isa DataFrame
    @test res.segments isa DataFrame
    @test nrow(res.segments) == 2
end

@testset "Somatic — Copy Number from Coverage" begin
    rng = MersenneTwister(19)
    cov = rand(rng, 100:200, 100)
    cn  = copy_number_from_coverage(cov; window=10)
    @test cn isa DataFrame
    @test hasproperty(cn, :log2_ratio)
    @test hasproperty(cn, :cn_call)
    @test nrow(cn) == 100
    @test all(c -> c in ("gain","loss","neutral"), cn.cn_call)

    # With reference
    ref_cov = rand(rng, 90:210, 100)
    cn_norm = copy_number_from_coverage(cov; window=10, reference_coverage=ref_cov)
    @test nrow(cn_norm) == 100
end

@testset "Somatic — SV Fusion Candidates" begin
    sv_df = DataFrame(
        chrom1 = ["chr1","chr2"],
        pos1   = [1000, 5000000],
        chrom2 = ["chr1","chr3"],
        pos2   = [5000, 8000000],
        sv_type = ["DEL","TRA"],
    )
    genes = DataFrame(
        chrom      = ["chr1","chr1","chr2","chr3"],
        gene_start = [500, 4500, 4999900, 7999900],
        gene_end   = [2000, 6000, 5000100, 8000100],
        gene_name  = ["GENE_A","GENE_B","GENE_C","GENE_D"],
    )
    fusions = sv_fusion_candidates(sv_df, genes; max_distance=200)
    @test fusions isa DataFrame
    @test nrow(fusions) >= 1
    @test hasproperty(fusions, :gene5p)
end

@testset "Somatic — Driver Enrichment" begin
    variants = DataFrame(
        gene        = ["TP53","TP53","KRAS","TTN","MYC","TP53"],
        consequence = ["missense_variant","stop_gained","missense_variant",
                       "synonymous_variant","missense_variant","frameshift_variant"],
        sample      = ["S1","S2","S1","S3","S1","S3"],
    )
    drivers = ["TP53","KRAS","MYC"]
    res = driver_enrichment(variants, drivers)
    @test res isa DataFrame
    @test hasproperty(res, :gene)
    @test hasproperty(res, :odds_ratio)
    @test hasproperty(res, :pvalue)
    @test nrow(res) > 0
end

@testset "Somatic — Variant Tier Classification" begin
    variants = DataFrame(
        consequence = ["stop_gained","missense_variant","synonymous_variant",
                       "frameshift_variant","intron_variant"],
    )
    tiered = variant_tier_classification(variants)
    @test tiered isa DataFrame
    @test hasproperty(tiered, :tier)
    @test hasproperty(tiered, :tier_label)
    @test tiered.tier[findfirst(==("stop_gained"), tiered.consequence)] == 1
    @test tiered.tier[findfirst(==("synonymous_variant"), tiered.consequence)] == 4
end

@testset "Somatic — Hotspot Scan" begin
    variants = DataFrame(
        gene     = ["TP53","TP53","TP53","KRAS","KRAS"],
        position = [248, 248, 248, 12, 12],
        sample   = ["S1","S2","S3","S1","S2"],
    )
    hs = somatic_hotspot_scan(variants; min_recurrence=2)
    @test hs isa DataFrame
    @test nrow(hs) >= 1
    @test all(hs.n_samples .>= 2)
end

@testset "Somatic — Mutational Spectrum 96" begin
    variants = DataFrame(
        trinucleotide_context = ["ACA","ACA","TCT","GCG"],
        ref = ["C","C","C","C"],
        alt = ["T","G","A","T"],
    )
    spec = mutational_spectrum_96(variants)
    @test spec isa DataFrame
    @test hasproperty(spec, :count)
    @test hasproperty(spec, :fraction)
    @test sum(spec.fraction) ≈ 1.0 atol=1e-8
end

@testset "CRISPR — Guide Design Heuristics" begin
    sequence = BioToolkit.DNASeq("ACGTACGTACGTACGTACGTAGG")
    guides = BioToolkit.design_guides(sequence)

    @test nrow(guides) >= 1
    @test guides.off_target_count[1] == -1
    @test ismissing(guides.specificity_score[1])

    filtered = BioToolkit.filter_guides(guides; min_specificity=0.9)
    @test nrow(filtered) == nrow(guides)

    ranked = BioToolkit.rank_guides(guides)
    @test hasproperty(ranked, :composite_score)
    @test all(isfinite, ranked.composite_score)
end

# =============================================================================
# IMMUNOLOGY  (immunology.jl new features)
# =============================================================================
@testset "Immunology — Somatic Hypermutation Rate" begin
    obs = ["ACDEFGCASSIREQFF", "ACDEGGCASSAREQFF", "ACDEFGCASKIRGQFF"]
    gl  = ["ACDEFGCASSIREQFF", "ACDEFGCASSIREQFF", "ACDEFGCASSIREQFF"]
    res = somatic_hypermutation_rate(obs, gl)
    @test res isa DataFrame
    @test hasproperty(res, :mutation_rate)
    @test nrow(res) == 3
    @test res.mutation_rate[1] ≈ 0.0 atol=1e-10   # identical to germline

    # With clone IDs
    clone_res = somatic_hypermutation_rate(obs, gl; clone_ids=["CL1","CL1","CL2"])
    @test clone_res isa NamedTuple
    @test clone_res.cells isa DataFrame
    @test clone_res.clones isa DataFrame
    @test nrow(clone_res.clones) == 2
end

@testset "Immunology — BCR Affinity Maturation Score" begin
    seqs = ["CASSIREQYF", "WWGACDEAGF", "MRKNASIREQYF"]
    res  = bcr_affinity_maturation_score(seqs)
    @test res isa DataFrame
    @test hasproperty(res, :hotspot_density)
    @test hasproperty(res, :maturation_score)
    @test nrow(res) == 3
    @test all(0.0 .<= res.maturation_score .<= 1.0)
end

@testset "Immunology — Antigen Specificity Clustering" begin
    cdr3s = ["CASSIREQYF","CASSIGEQYF","CASSIREQYF","LARGHDQETQYF","LARGADQETQYF"]
    res   = antigen_specificity_clustering(cdr3s; similarity_threshold=2.0, min_cluster_size=2)
    @test res isa NamedTuple
    @test res.cells isa DataFrame
    @test res.cluster_summary isa DataFrame
    @test hasproperty(res.cells, :cluster_id)
    # The first 3 identical/near-identical sequences should cluster together
    ids = res.cells.cluster_id[1:3]
    @test all(==(ids[1]), ids)
end

@testset "Immunology — Clonotype Trajectory" begin
    rng  = MersenneTwister(23)
    clones = repeat(["CL_A","CL_B","CL_C"], inner=10)
    tps    = repeat(["T1","T2","T3","T4","T5"], outer=6)

    # Make CL_A expand over time
    traj = clonotype_trajectory(clones, tps; min_max_frequency=0.0)
    @test traj isa DataFrame
    @test hasproperty(traj, :delta_first_last)
    @test hasproperty(traj, :expansion_class)
    @test all(c -> c in ("expanding","contracting","stable"), traj.expansion_class)
end

@testset "Immunology — V-Gene Family Usage" begin
    vgenes = ["TRBV12-3","TRBV12-4","TRBV20-1","TRBV12-1","IGHV1-2","TRBV20-1"]
    res    = v_gene_family_usage(vgenes)
    @test res isa DataFrame
    @test hasproperty(res, :v_family)
    @test hasproperty(res, :frequency)
    @test sum(res.frequency) ≈ 1.0 atol=1e-8

    # With groups
    groups = ["A","A","A","B","B","B"]
    res_g  = v_gene_family_usage(vgenes; group_col=groups)
    @test hasproperty(res_g, :group)
end

@testset "Immunology — HLA Type Enrichment" begin
    rng = MersenneTwister(29)
    n = 60
    clonotypes = [rand(rng, ["CL_A","CL_B","CL_C","CL_D"]) for _ in 1:n]
    hla_alleles = [rand(rng, ["HLA-A*02:01","HLA-A*01:01","HLA-B*07:02"]) for _ in 1:n]
    res = hla_type_enrichment(clonotypes, hla_alleles; top_fraction=0.25)
    @test res isa DataFrame
    @test hasproperty(res, :hla_allele)
    @test hasproperty(res, :pvalue)
    @test all(0.0 .<= res.padj .<= 1.0)
end

@testset "Immunology — CDR3 Physicochemical Properties" begin
    seqs = [BioToolkit.AASeq("CASSIREQYF"), BioToolkit.AASeq("CARTGDKTEV"), BioToolkit.AASeq("CAASMDSNYQLIW")]
    res  = cdr3_physicochemical_properties(seqs)
    @test res isa DataFrame
    @test hasproperty(res, :net_charge)
    @test hasproperty(res, :hydrophobicity)
    @test hasproperty(res, :aromaticity)
    @test hasproperty(res, :instability_idx)
    @test nrow(res) == 3
    @test all(isfinite, res.hydrophobicity)
end

@testset "Immunology — Network Centrality" begin
    cdr3s = ["CASSIREQYF","CASSIGEQYF","CASSGGEQYF","LARGHDQETQYF","CARGRDQETQYF"]
    res   = network_centrality_clonotypes(cdr3s; similarity_threshold=3.0)
    @test res isa DataFrame
    @test hasproperty(res, :degree)
    @test hasproperty(res, :hub_score)
    @test nrow(res) == 5
    @test all(res.degree .>= 0)
end

@testset "Immunology — Convergent Clonotype Detection" begin
    seqs  = ["CASSIREQYF","CASSIREQYF","CASSIGEQYF","LARGVDQETQYF","LARGVDQETQYF"]
    vgenes = ["TRBV12","TRBV14","TRBV12","TRBV20","TRBV22"]
    jgenes = ["TRBJ1","TRBJ1","TRBJ2","TRBJ2","TRBJ3"]
    res   = convergent_clonotype_detection(seqs, vgenes, jgenes; max_hamming=1, min_group_size=2)
    @test res isa DataFrame
    @test hasproperty(res, :is_convergent)
    @test nrow(res) >= 1
end

@testset "Immunology — Lineage Tree from Clones" begin
    germline = BioToolkit.AASeq("CASSIREQYF")
    seqs     = [BioToolkit.AASeq("CASSIREQYF"), BioToolkit.AASeq("CASSIGEQYF"), BioToolkit.AASeq("CASSIEEQYF"), BioToolkit.AASeq("CASTIGAQYF")]
    tree     = lineage_tree_from_clones(seqs, germline)
    @test tree isa NamedTuple
    @test tree.edges isa DataFrame
    @test hasproperty(tree.edges, :parent)
    @test hasproperty(tree.edges, :child)
    @test hasproperty(tree.edges, :hamming_dist)
    @test nrow(tree.edges) == 4   # n_seqs edges (germline + 4 seqs → 4 MST edges)
    @test isfinite(tree.mst_distance_sum)
end

@testset "Immunology — Epitope Binding Profile" begin
    peptides = [BioToolkit.AASeq("GILGFVFTL"), BioToolkit.AASeq("NLVPMVATV"), BioToolkit.AASeq("YVNVNMGLK")]
    known    = [BioToolkit.AASeq("GILGFVFTL"), BioToolkit.AASeq("NLVPMVATV")]
    res      = epitope_binding_profile(peptides, known; allele="HLA-A*02:01", max_hamming=0)
    @test res isa DataFrame
    @test hasproperty(res, :combined_binding_score)
    @test hasproperty(res, :nearest_epitope)
    @test nrow(res) == 3
    @test res.predicted_binder[1] == true  # exact match to known epitope
end

# =============================================================================
# DEEP LEARNING  (deeplearning.jl new features)
# =============================================================================
@testset "DeepLearning — Cell Type Denoising (MAGIC-like)" begin
    rng   = MersenneTwister(31)
    X     = abs.(rand(rng, Int, 20, 50) .% 100) .+ 1
    res   = cell_type_denoising(X; k=5, t=2, n_pcs=5, backend=:cpu)
    @test res isa NamedTuple
    @test res.denoised isa Matrix{Float64}
    @test size(res.denoised) == (20, 50)
    @test all(isfinite, res.denoised)
    @test res.diffusion_operator isa Matrix{Float64}
    @test size(res.diffusion_operator) == (50, 50)
end

@testset "DeepLearning — Sparse Autoencoder Features" begin
    rng = MersenneTwister(37)
    X   = abs.(rand(rng, Int, 30, 40) .% 50) .+ 1
    res = sparse_autoencoder_features(X; n_features=16, sparsity_k=4, n_iter=50, seed=1, backend=:cpu)
    @test res isa NamedTuple
    @test res.features isa Matrix{Float64}
    @test size(res.features) == (30, 16)
    @test res.feature_activations isa Matrix{Float64}
    @test size(res.feature_activations) == (40, 16)
    @test isfinite(res.reconstruction_loss)
    # Sparsity: each row should have at most sparsity_k non-zeros
    for row in eachrow(res.feature_activations)
        @test count(!iszero, row) <= 4
    end
end

@testset "DeepLearning — Trajectory Neural ODE" begin
    rng = MersenneTwister(41)
    S   = abs.(rand(rng, Int, 15, 30) .% 80) .+ 5
    U   = abs.(rand(rng, Int, 15, 30) .% 40) .+ 2
    res = trajectory_neural_ode(S, U; n_latent=5, n_steps=10, dt=0.05, backend=:cpu)
    @test res isa NamedTuple
    @test res.pseudotime isa Vector{Float64}
    @test length(res.pseudotime) == 30
    @test all(0.0 .<= res.pseudotime .<= 1.0)
    @test res.latent_velocity isa Matrix{Float64}
    @test size(res.trajectory, 3) == 11    # n_steps + 1
end

@testset "DeepLearning — Multimodal WNN Embedding" begin
    rng = MersenneTwister(43)
    n_cells = 25
    rna  = abs.(rand(rng, Int, 20, n_cells) .% 100) .+ 1
    atac = abs.(rand(rng, Int, 10, n_cells) .% 50)  .+ 1
    res  = multimodal_wnn_embedding([rna, atac]; n_latent=8, n_neighbors=5, seed=1, backend=:cpu)
    @test res isa NamedTuple
    @test res.latent isa Matrix{Float64}
    @test size(res.latent, 1) == n_cells
    @test res.modality_weights isa Matrix{Float64}
    @test size(res.modality_weights) == (n_cells, 2)
    @test all(isfinite, res.latent)
    # Weights should sum to ~1 per cell
    for row in eachrow(res.modality_weights)
        @test sum(row) ≈ 1.0 atol=1e-6
    end
end

@testset "DeepLearning — Protein Sequence Embedding" begin
    seqs = BioToolkit.AASeq.( ["MKTLLILAVLCL", "ACDEFGHIKLMNPQRSTVWY", "CASSIREQYFP"] )
    emb  = protein_sequence_embedding(seqs; dim=32, k=3, include_properties=true)
    @test emb isa Matrix{Float64}
    @test size(emb) == (3, 32)
    @test all(isfinite, emb)
    # Each row should be normalised
    for row in eachrow(emb)
        @test norm(row) ≈ 1.0 atol=1e-6
    end
end

@testset "DeepLearning — Zero-Shot Cell Annotation" begin
    rng = MersenneTwister(47)
    X   = abs.(rand(rng, Int, 10, 20) .% 100) .+ 1
    genes = ["G$i" for i in 1:10]
    archetypes = Dict(
        "T_cell"  => ["G1","G2","G3"],
        "B_cell"  => ["G4","G5","G6"],
        "NK_cell" => ["G7","G8","G9"],
    )
    res = zero_shot_cell_annotation(X, genes, archetypes)
    @test res isa NamedTuple
    @test res.cells isa DataFrame
    @test hasproperty(res.cells, :predicted_label)
    @test nrow(res.cells) == 20
    @test all(l -> l in ["T_cell","B_cell","NK_cell"], res.cells.predicted_label)
    @test size(res.score_matrix) == (20, 3)
end

@testset "DeepLearning — GNN Gene Regulatory Network" begin
    rng = MersenneTwister(53)
    X   = abs.(rand(rng, Int, 20, 40) .% 80) .+ 2
    res = gene_regulatory_network_gnn(X; k_neighbors=5, n_propagation=2, top_k_edges=30, n_pcs=5, backend=:cpu)
    @test res isa NamedTuple
    @test res.edges isa DataFrame
    @test hasproperty(res.edges, :source_gene)
    @test hasproperty(res.edges, :propagated_weight)
    @test nrow(res.edges) <= 30
end

@testset "DeepLearning — Self-Supervised Pretraining" begin
    rng = MersenneTwister(59)
    X   = abs.(rand(rng, Int, 25, 30) .% 60) .+ 3
    res = self_supervised_pretraining(X; mask_fraction=0.15, n_latent=8, n_iter=30, seed=1, backend=:cpu)
    @test res isa NamedTuple
    @test res.latent isa Matrix{Float64}
    @test size(res.latent, 1) == 30
    @test size(res.latent, 2) == 8
    @test !isempty(res.reconstruction_loss_history)
end

@testset "DeepLearning — Cell Cycle Regression" begin
    rng = MersenneTwister(61)
    X   = abs.(rand(rng, Int, 15, 20) .% 50) .+ 1
    genes = ["G$i" for i in 1:15]
    s_genes   = ["G1","G2","G3"]
    g2m_genes = ["G4","G5","G6"]
    res = cell_cycle_regression(X, genes; s_genes=s_genes, g2m_genes=g2m_genes)
    @test res isa NamedTuple
    @test res.corrected_counts isa Matrix{Float64}
    @test size(res.corrected_counts) == (15, 20)
    @test res.s_score isa Vector{Float64}
    @test length(res.g2m_score) == 20
    @test all(p -> p in ["S","G2M","G1"], res.phase)
end

@testset "DeepLearning — Deep Factorization (MOFA+)" begin
    rng = MersenneTwister(67)
    n_cells = 20
    rna  = abs.(rand(rng, Int, 15, n_cells) .% 80) .+ 1
    atac = abs.(rand(rng, Int, 10, n_cells) .% 40) .+ 1
    res  = deep_factorization_embedding([rna, atac]; n_factors=5, n_iter=50, seed=1, backend=:cpu)
    @test res isa NamedTuple
    @test res.factors isa Matrix{Float64}
    @test size(res.factors, 1) == n_cells
    @test length(res.loadings) == 2
    @test size(res.loadings[1]) == (15, 5)
    @test isfinite(res.reconstruction_error)
end

# =============================================================================
# COEVOLUTION  (coevolution.jl new features)
# =============================================================================
@testset "Coevolution — Mutual Information Contacts" begin
    msa = _make_msa(50, 12)
    cm  = mutual_information_contacts(msa; min_separation=2, apc=true)
    @test cm isa ContactMap
    @test size(cm.scores, 1) == size(cm.scores, 2)
    @test cm.scores ≈ cm.scores'    # symmetric
    @test all(cm.scores[i,i] == 0 for i in axes(cm.scores,1))
    @test all(cm.scores .>= 0)
end

@testset "Coevolution — Direct Information Contacts" begin
    msa = _make_msa(50, 12)
    cm  = direct_information_contacts(msa; min_separation=2)
    @test cm isa ContactMap
    @test size(cm.scores, 1) == size(cm.scores, 2)
    @test cm.scores ≈ cm.scores'
    @test maximum(cm.scores) <= 1.0 + 1e-8
end

@testset "Coevolution — Column Conservation Scores" begin
    msa  = _make_msa(30, 10)
    cons = column_conservation_scores(msa)
    @test cons isa Vector{Float64}
    @test length(cons) == 10
    @test all(0.0 .<= cons .<= 1.0)
    # Identical column should have conservation = 1
    identical_seqs = MultipleSequenceAlignment([
        "ACGT", "ACGT", "ACGT", "ACGT"
    ])
    cons_identical = column_conservation_scores(identical_seqs)
    @test all(cons_identical .≈ 1.0)
end

@testset "Coevolution — Sequence Logo Entropy" begin
    msa  = _make_msa(20, 8)
    logo = sequence_logo_entropy(msa)
    @test logo isa NamedTuple
    @test hasproperty(logo, :position)
    @test hasproperty(logo, :entropy)
    @test hasproperty(logo, :information_content)
    @test length(logo.position) == 8
    @test all(logo.entropy .>= 0)
end

@testset "Coevolution — Evolutionary Coupling Network" begin
    msa = _make_msa(60, 12)
    model = fit_pseudolikelihood_model(msa; regularization=0.05)
    scores_mat = compute_contact_scores(model; apc=true, min_separation=1)
    cm = ContactMap(scores_mat, collect(1:12))
    # Normalise for network construction
    cm_norm = ContactMap(
        let s = cm.scores; max_v = maximum(s); max_v > 0 ? s ./ max_v : s end,
        cm.residue_ids
    )
    net = evolutionary_coupling_network(cm_norm; score_threshold=0.0, min_separation=1)
    @test net isa NamedTuple
    @test net.edges isa NamedTuple
    @test length(net.degree) == 12
    @test net.n_edges >= 0
end

@testset "Coevolution — Contact Enrichment Statistics" begin
    l = 10
    scores = zeros(Float64, l, l)
    for i in 1:l-1, j in (i+1):l
        if abs(i-j) >= 2
            scores[i,j] = scores[j,i] = rand() * 0.8 + 0.1
        end
    end
    cm = ContactMap(scores, collect(1:l))
    true_contacts = zeros(Float64, l, l)
    true_contacts[1,6] = true_contacts[6,1] = 1.0
    true_contacts[2,7] = true_contacts[7,2] = 1.0
    stats = contact_enrichment_statistics(cm, true_contacts; top_fractions=[0.5,1.0], min_separation=2)
    @test stats isa NamedTuple
    @test haskey(stats.precision_at_k, 0.5)
    @test haskey(stats.precision_at_k, 1.0)
    @test isfinite(stats.ppv_auc)
    @test isfinite(stats.mcc)
end

@testset "Coevolution — Gap Analysis" begin
    msa = MultipleSequenceAlignment([
        "AC-DEFG",
        "ACDDEFG",
        "A--DEFG",
        "ACDDE-G",
    ])
    ga = gap_analysis(msa)
    @test ga isa NamedTuple
    @test length(ga.column_gap_fraction) == 7
    @test ga.column_gap_fraction[3] ≈ 0.5 atol=0.01   # col 3: '-' in 2/4 seqs
    @test length(ga.sequence_gap_fraction) == 4
    @test 0.0 <= ga.total_gap_fraction <= 1.0
    @test length(ga.gap_blocks) == 4
end

@testset "Coevolution — Contact Precision-Recall" begin
    l = 10
    scores = zeros(Float64, l, l)
    for i in 1:l-1, j in (i+1):l
        abs(i-j) >= 3 && (scores[i,j] = scores[j,i] = rand())
    end
    cm = ContactMap(scores, collect(1:l))
    true_contacts = zeros(Float64, l, l)
    true_contacts[1,8] = true_contacts[8,1] = 1.0
    pr = contact_precision_recall(cm, true_contacts; min_separation=3, n_points=10)
    @test pr isa NamedTuple
    @test length(pr.thresholds) == 10
    @test length(pr.precision) == 10
    @test length(pr.recall) == 10
    @test 0.0 <= pr.auc_pr <= 1.0
end

@testset "Coevolution — Phylogenetic Correction (Henikoff)" begin
    msa = _make_msa(20, 8)
    res = phylogenetic_correction(msa; method=:henikoff)
    @test res isa NamedTuple
    @test length(res.weights) == 20
    @test all(res.weights .> 0)
    @test isfinite(res.effective_sequences)
    # Clustering method
    res2 = phylogenetic_correction(msa; method=:clustering, identity_threshold=0.9)
    @test length(res2.weights) == 20
end

@testset "Coevolution — Positional Covariation Matrix" begin
    msa = _make_msa(40, 10)
    res = positional_covariation_matrix(msa; metric=:pearson)
    @test res isa NamedTuple
    @test res.covariation isa Matrix{Float64}
    @test size(res.covariation) == (10, 10)
    @test res.covariation ≈ res.covariation'   # symmetric
    @test all(abs.(diag(res.covariation)) .<= 1.0 + 1e-6)   # |correlation| ≤ 1
end

@testset "Coevolution — Alignment Quality Report" begin
    msa = _make_msa(30, 12)
    rep = alignment_quality_report(msa)
    @test rep isa NamedTuple
    @test rep.n_sequences == 30
    @test rep.alignment_length == 12
    @test rep.effective_sequences > 0
    @test 0.0 <= rep.mean_pairwise_identity <= 1.0
    @test 0.0 <= rep.mean_gap_fraction <= 1.0
    @test 0.0 <= rep.mean_conservation <= 1.0
    @test length(rep.conservation_scores) == 12
    @test length(rep.entropy_per_column) == 12
end
