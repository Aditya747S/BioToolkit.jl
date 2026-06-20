using Test
using Random
using Statistics
using DataFrames
using BioToolkit

@testset "Domain feature expansion" begin
    @testset "Sequence utilities" begin
        dna = BioToolkit.DNASeq("ACGTACGT")

        oligos = BioToolkit.oligonucleotide_frequency(dna; k=2)
        @test oligos["AC"] == 2
        @test oligos["CG"] == 2
        @test oligos["GT"] == 2
        @test oligos["TA"] == 1

        normalized = BioToolkit.oligonucleotide_frequency(dna; k=2, normalize=true)
        @test isapprox(sum(values(normalized)), 1.0; atol=1e-12)

        @test isapprox(BioToolkit.sequence_entropy(dna), 2.0; atol=1e-12)

        windowed = BioToolkit.sliding_window_gc_content(dna; window=4, step=2)
        @test windowed.start_positions == [1, 3, 5]
        @test windowed.stop_positions == [4, 6, 8]
        @test all(isapprox.(windowed.gc_content, 0.5; atol=1e-12))
    end

    @testset "Spatial advanced" begin
        counts = [
            10 9 1 0;
            9 8 1 0;
            1 1 9 10;
            0 1 8 9;
            5 4 5 4
        ]
        gene_ids = ["L1", "R1", "L2", "R2", "HK"]
        cell_ids = ["c1", "c2", "c3", "c4"]
        coords = [0.0 0.0; 0.1 0.0; 2.0 2.0; 2.1 2.0]
        sce = BioToolkit.SingleCellExperiment(counts, gene_ids, cell_ids; spatial_coords=coords)

        domains = BioToolkit.bayesspace_like_domains(sce; n_domains=2, n_pcs=3)
        @test nrow(domains.assignments) == 4

        spagc = BioToolkit.spagc_like_domains(sce; n_domains=2, n_pcs=3, n_iter=3, k=2)
        @test nrow(spagc) == 4

        svg = BioToolkit.sparkx_spatial_de(sce; k=2, n_permutations=20, seed=4)
        @test nrow(svg) == length(gene_ids)

        sde = BioToolkit.spatialde_gp_de(sce; length_scale=0.8, k=2)
        @test nrow(sde) == length(gene_ids)

        traj = BioToolkit.spatial_trajectory_graph(sce; n_components=3, k=2)
        @test length(traj.pseudotime) == length(cell_ids)

        labels = [1, 1, 2, 2]
        niche = BioToolkit.niche_weighted_communication(sce, labels)
        @test "weighted_score" in names(niche)

        refined = BioToolkit.spatial_markov_refine_domains(sce, labels; n_iter=4, k=2, threaded=true)
        @test nrow(refined) == length(cell_ids)

        perm = BioToolkit.spatial_lr_permutation_test(sce, labels; n_permutations=8, seed=9, threaded=true)
        @test "padj" in names(perm)

        seg = BioToolkit.cell2location_like_segmentation(rand(3, 5), rand(2, 5); cell_type_names=["T", "B"])
        @test nrow(seg) == 3
        @test "dominant_celltype" in names(seg)
    end

    @testset "Immunology" begin
        seq = BioToolkit.AASeq("MRCASSPGGGGFGK")
        @test !ismissing(BioToolkit.extract_cdr3(seq))

        clono = BioToolkit.clonotype_table([BioToolkit.AASeq("CASSIRSSYEQYF"), BioToolkit.AASeq("CASSIRSSYEQYF"), BioToolkit.AASeq("CASSLGQNTLYF")])
        @test nrow(clono) >= 1

        vdb = Dict("TRBV1" => "CASSIR", "TRBV2" => "LGQNTL")
        jdb = Dict("TRBJ1" => "EQYF", "TRBJ2" => "TLYF")
        vdj = BioToolkit.assign_vdj_segments([BioToolkit.AASeq("CASSIRSSYEQYF"), BioToolkit.AASeq("CASSLGQNTLYF")], vdb, jdb)
        @test nrow(vdj) == 2

        iso = BioToolkit.isotype_switching_summary(["c1", "c1", "c2"], ["IgM", "IgG1", "IgM"])
        @test "switches" in names(iso)

        usage = BioToolkit.germline_usage_bias(["TRBV1", "TRBV1", "TRBV2"])
        @test nrow(usage) == 2

        b = BioToolkit.bepipred_like_scores(BioToolkit.AASeq("ACDEFGHIKLMNPQRSTVWY"))
        @test nrow(b) == 20

        mhc = BioToolkit.mhcflurry_like_scores([BioToolkit.AASeq("SLYNTVATL"), BioToolkit.AASeq("AAAAAAA")])
        @test nrow(mhc) == 2

        bcr = BioToolkit.bcr_affinity_maturation_score([BioToolkit.DNASeq("ATGCGTATGC"), BioToolkit.DNASeq("AATTGGCCAA")])
        @test nrow(bcr) == 2

        div = BioToolkit.repertoire_diversity_metrics(["c1", "c1", "c2", "c3"])
        @test div.richness == 3

        expn = BioToolkit.clonal_expansion_test(["a", "a", "b", "c", "a"], ["case", "case", "ctrl", "ctrl", "ctrl"]; case_group="case")
        @test nrow(expn) == 3

        len_spec = BioToolkit.cdr3_length_spectrum([BioToolkit.AASeq("CASSIRSSYEQYF"), BioToolkit.AASeq("MRCASSPGGGGFGK"), BioToolkit.AASeq("NOHITSEQ")])
        @test nrow(len_spec) >= 1

        dist = BioToolkit.tcrdist_like_matrix([BioToolkit.AASeq("CASSIRSSYEQYF"), BioToolkit.AASeq("CASSLGQNTLYF"), BioToolkit.AASeq("CASSPGQETQYF")])
        @test size(dist) == (3, 3)

        igb = BioToolkit.igblast_assign("vdj.fa"; dry_run=true)
        @test igb.status == :planned

        manifest = BioToolkit.imgt_highvquest_manifest(Dict("s1" => "a.fa", "s2" => "b.fa"))
        @test nrow(manifest) == 2

        paired = BioToolkit.paired_clonotype_table([BioToolkit.AASeq("CASSA"), BioToolkit.AASeq("CASSA"), BioToolkit.AASeq("CASSB")], [BioToolkit.AASeq("CSARQ"), BioToolkit.AASeq("CSARQ"), BioToolkit.AASeq("CSART")])
        @test nrow(paired) >= 1

        dvj = BioToolkit.differential_vj_usage(["TRBV1", "TRBV1", "TRBV2", "TRBV2"], ["TRBJ1", "TRBJ1", "TRBJ2", "TRBJ2"], ["case", "ctrl", "case", "ctrl"]; case_group="case")
        @test nrow(dvj) == 2
    end

    @testset "Multi-omics advanced" begin
        Random.seed!(11)
        assay1 = randn(12, 20)
        assay2 = randn(12, 15)
        mofa = BioToolkit.mofa_plus_integration([assay1, assay2]; n_factors=4)
        @test size(mofa.factors, 2) == 4

        cca = BioToolkit.cca_integration(randn(30, 8), randn(30, 6); n_components=3)
        @test size(cca.x_scores, 2) == 3

        anchors = BioToolkit.anchor_based_integration(randn(25, 10), randn(12, 10); n_components=5, k=2)
        @test nrow(anchors) == 24

        totalvi = BioToolkit.totalvi_like_integration(rand(40, 18), rand(8, 18); n_latent=6)
        @test size(totalvi.latent, 2) == 6

        mofa_em = BioToolkit.mofa_plus_em([randn(10, 12), randn(8, 12)]; n_factors=3, n_iter=20, backend=:cpu, threaded=true)
        @test size(mofa_em.factors, 2) == 3

        scca = BioToolkit.sparse_cca_integration(randn(40, 10), randn(40, 7); n_components=3, sparsity=0.3, threaded=true)
        @test size(scca.x_scores, 2) == 3

        cdf = DataFrame(treat=[1, 1, 0, 0, 1, 0], outcome=[4.2, 4.5, 3.1, 3.0, 4.0, 3.3], age=[50, 52, 48, 47, 51, 49])
        ate = BioToolkit.causal_ate_regression(cdf, :treat, :outcome; covariates=[:age])
        @test isfinite(ate.ate)
    end

    @testset "Long-read additions" begin
        cons = BioToolkit.isoseq_consensus(["ACGTAC", "ACGTTC", "TTTTAA", "TTTTAT"]; n_clusters=2)
        @test nrow(cons) >= 1

        mm2 = BioToolkit.minimap2_align("ref.fa", "reads.fq")
        @test mm2.status == :planned

        ga = BioToolkit.graphaligner_align("graph.gfa", "reads.fq")
        @test ga.status == :planned

        snf = BioToolkit.sniffles_call("reads.bam")
        @test snf.status == :planned

        svim = BioToolkit.svim_call("reads.bam", "ref.fa")
        @test svim.status == :planned

        wh = BioToolkit.whatshap_phase("calls.vcf", "reads.bam")
        @test wh.status == :planned

        phase = BioToolkit.phase_reads_by_alleles([0 1 -1; 0 1 1; 1 0 0; 1 0 -1])
        @test length(phase.read_haplotype) == 4

        olc = BioToolkit.overlap_layout_consensus(["ACGTACGT", "ACGTGG", "GGTTTT"]; min_overlap=2)
        @test ncodeunits(olc.consensus) >= 6

        n50 = BioToolkit.read_n50(["AAAA", "AAAAAAAA", "AA"])
        @test n50 == 8

        og = BioToolkit.overlap_graph_table(["ACGTACGT", "ACGTGG", "GGTTTT"]; min_overlap=2, threaded=true)
        @test nrow(og) >= 1

        segs = DataFrame(read_id=["r1", "r1", "r2", "r2"], chrom=["chr1", "chr1", "chr1", "chr2"], start=[100, 500, 200, 300], stop=[150, 550, 250, 350])
        pc = BioToolkit.porec_contact_table(segs)
        @test nrow(pc) >= 1

        cm = BioToolkit.omnic_contact_matrix(pc; chrom="chr1", bin_size=100)
        @test size(cm.matrix, 1) == size(cm.matrix, 2)
    end

    @testset "Deep-learning style utilities" begin
        emb = BioToolkit.scvi_like_embedding(rand(60, 20); n_latent=5, backend=:cpu)
        @test size(emb.latent, 2) == 5

        expr = rand(5, 10)
        gene_ids = ["G1", "G2", "G3", "G4", "G5"]
        markers = Dict("T" => ["G1", "G2"], "B" => ["G3", "G4"])
        ca = BioToolkit.cellassign_like_mapping(expr, gene_ids, markers)
        @test length(ca.predicted_label) == 10

        gf = BioToolkit.geneformer_like_embedding(BioToolkit.DNASeq.( ["ACGTACGT", "TTTTGGGG"] ); k=3, dim=16, threaded=true)
        @test size(gf) == (2, 16)

        gpt = BioToolkit.scgpt_like_embedding(BioToolkit.DNASeq.( ["ACGTACGT", "TTTTGGGG"] ); token_dim=24, threaded=true)
        @test size(gpt) == (2, 24)

        bert = BioToolkit.scbert_like_embedding(BioToolkit.DNASeq.( ["ACGTACGT", "TTTTGGGG"] ); dim=20, max_len=32, threaded=true)
        @test size(bert) == (2, 20)

        grn = BioToolkit.attention_grn(rand(20, 15); top_k=10)
        @test nrow(grn) > 0

        bc = BioToolkit.batch_corrected_latent(rand(50, 12), ["b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2", "b3", "b3", "b3", "b3"]; n_latent=4, backend=:cpu)
        @test size(bc.latent, 2) == 4

        con = BioToolkit.contrastive_cell_embedding(rand(30, 14); n_latent=5, dropout=0.2, seed=2, backend=:cpu)
        @test size(con.latent, 2) == 5

        sg = BioToolkit.scgen_like_perturbation(rand(20, 8), rand(20, 8), rand(20, 5); n_latent=4, backend=:cpu)
        @test size(sg.predicted_treated, 2) == 5

        A = rand(10, 10)
        A = 0.5 .* (A .+ A')
        gl = BioToolkit.graphsca_label_transfer(rand(15, 10), A, ["T", "T", "unknown", "B", "B", "unknown", "T", "B", "unknown", "T"])
        @test length(gl.predicted_label) == 10

        if Base.find_package("Flux") !== nothing
            import Flux
            ae = BioToolkit.flux_autoencoder_embedding(rand(40, 12); latent_dim=4, hidden_dim=8, epochs=1, backend=:cpu)
            @test size(ae.latent, 2) == 4

            clf = BioToolkit.flux_mlp_classifier(rand(12, 20), vcat(fill("A", 10), fill("B", 10)); hidden_dim=8, epochs=1, backend=:cpu)
            @test length(clf.predicted_label) == 20
        else
            @test_throws ArgumentError BioToolkit.flux_autoencoder_embedding(rand(20, 10); epochs=1)
        end
    end

    @testset "Epigenomics advanced" begin
        states = BioToolkit.chromhmm_like_segmentation(randn(100, 4); n_states=5)
        @test length(states.state) == 100

        mat = abs.(randn(20, 20))
        hic = BioToolkit.hic_ice_normalize(mat; max_iter=15)
        @test size(hic.normalized) == (20, 20)

        kr = BioToolkit.hic_kr_normalize(mat; max_iter=25)
        @test size(kr.normalized) == (20, 20)

        ab = BioToolkit.hic_ab_compartments(mat)
        @test length(ab.compartment) == 20

        fp = BioToolkit.tobias_like_footprints(rand(200), [40, 80, 120])
        @test nrow(fp) == 3

        bismark = BioToolkit.parse_bismark_coverage([
            "chr1\t100\t100\t75.0\t3\t1",
            "chr1\t101\t101\t20.0\t1\t4"
        ])
        @test nrow(bismark.calls) == 2
        @test nrow(bismark.summary) >= 1

        ins = BioToolkit.insulation_score(abs.(randn(24, 24)); window=4)
        @test nrow(ins) == 24

        bnd = BioToolkit.compartment_boundary_candidates(ins; q=0.3)
        @test "is_boundary" in names(bnd)
    end

    @testset "Proteomics advanced" begin
        dia = BioToolkit.dia_like_quantification(rand(12, 8), [:A, :A, :A, :A, :B, :B, :B, :B])
        @test nrow(dia) == 12

        ph = BioToolkit.phosphosite_localization(BioToolkit.AASeq("AASTY"), Dict(3 => 0.8, 4 => 0.2))
        @test ph.confidence == "class-I"

        gly = BioToolkit.glycoproteomics_motif_table([BioToolkit.AASeq("NVTAA"), BioToolkit.AASeq("NPSTA")])
        @test nrow(gly) == 2

        embed = Dict("A" => [1.0, 0.0], "C" => [0.5, 0.5], "D" => [0.0, 1.0])
        proj = BioToolkit.project_peptides_to_structure([BioToolkit.AASeq("ACD"), BioToolkit.AASeq("DA")], embed)
        @test size(proj) == (2, 2)

        top3 = BioToolkit.protein_inference_top3(DataFrame(protein=["P1", "P1", "P1", "P2"], intensity=[10.0, 6.0, 4.0, 9.0]))
        @test nrow(top3) == 2

        ptm = BioToolkit.ptm_site_enrichment(["S10", "S10", "T20"], ["S10", "T20", "Y5", "S10", "S3"])
        @test nrow(ptm) >= 1
    end

    @testset "Clinical advanced" begin
        calls = DataFrame(gene=["CYP2D6", "CYP2C19"], variant=["rs1065852", "rs4244285"])
        pgx = BioToolkit.pharmacogenomics_star_alleles(calls)
        @test nrow(pgx) == 2

        person = DataFrame(person_id=[1, 2], age=[60, 45])
        visit = DataFrame(person_id=[1, 1, 2], visit_id=[1, 2, 3])
        cond = DataFrame(person_id=[1, 2, 2], condition_id=[10, 11, 12])
        omop = BioToolkit.omop_visit_summary(person, visit, cond)
        @test nrow(omop) == 2

        synth = BioToolkit.synthpop_like_cohort(DataFrame(a=[1.0, 2.0, 3.0], b=["x", "y", "z"]); n=5, seed=3)
        @test nrow(synth) == 5

        clinical = DataFrame(age=[50, 80], ecog=[0, 2], biomarker=[1.0, 0.0])
        score = BioToolkit.trial_suitability_scores(clinical; biomarker_cols=[:biomarker])
        @test nrow(score) == 2
        @test "eligible" in names(score)

        pheno = BioToolkit.cpic_metabolizer_phenotype(["*1/*4", "*1/*1"]; gene="CYP2D6")
        @test nrow(pheno) == 2

        rec = BioToolkit.pharmgkb_like_recommendations(DataFrame(sample_id=["s1"], gene=["CYP2D6"], diplotype=["*1/*4"]); sample_col=:sample_id)
        @test nrow(rec) == 1

        psm = BioToolkit.propensity_score_match(DataFrame(treatment=[1, 1, 0, 0, 0], age=[60.0, 55.0, 58.0, 40.0, 70.0], ecog=[1.0, 0.0, 1.0, 2.0, 0.0]); treatment_col=:treatment, covariates=[:age, :ecog], ratio=1)
        @test nrow(psm.matches) >= 1
    end

    @testset "Metagenomics advanced" begin
        bins = BioToolkit.mag_bin_contigs(rand(20, 6), rand(20); n_bins=4)
        @test nrow(bins) == 20

        db = Dict("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" => "Bacteria")
        tax = BioToolkit.kraken_like_classify([BioToolkit.DNASeq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), BioToolkit.DNASeq("TTTT")], db; k=31)
        @test nrow(tax) == 2

        viral = BioToolkit.viral_contig_scores([BioToolkit.DNASeq("ATATATAT"), BioToolkit.DNASeq("GCGCGCTATAAATAAA")])
        @test nrow(viral) == 2

        humann = BioToolkit.humann_like_pathways(DataFrame(feature=["U1", "U2"], abundance=[1.0, 2.0]), Dict("U1" => "P1", "U2" => "P1"))
        @test nrow(humann) == 1

        strain = BioToolkit.strainge_like_variants(DataFrame(ref_count=[20, 5], alt_count=[1, 10]); min_depth=5, min_vaf=0.1)
        @test nrow(strain) == 2
        @test "is_variant" in names(strain)

        lca = BioToolkit.lca_taxonomy_from_votes([["k__Bacteria;p__Firmicutes;c__Bacilli", "k__Bacteria;p__Firmicutes;c__Bacilli"], ["k__Bacteria;p__Proteobacteria", "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria"]])
        @test nrow(lca) == 2

        hap = BioToolkit.strain_haplotype_profile(DataFrame(sample=["s1", "s1", "s2"], position=[100, 105, 100], allele=["A", "G", "T"]))
        @test nrow(hap) == 2
    end

    @testset "Pipeline cloud and distributed" begin
        plan = BioToolkit.terra_workspace_plan(DataFrame(sample_id=["s1", "s2"]); workspace_name="ws")
        @test nrow(plan) == 2

        manifest = BioToolkit.anvil_workspace_manifest(["gs://bucket/a.bam", "gs://bucket/b.bam"])
        @test nrow(manifest) == 2

        mapped = BioToolkit.distributed_pipeline_map(x -> x^2, 1:8; nworkers=2)
        @test mapped == [x^2 for x in 1:8]

        cnode = BioToolkit.containerized_node(:noop, () -> 1; image="ghcr.io/biotoolkit/base:latest")
        @test haskey(cnode.inputs, :container_image)

        pg = BioToolkit.PipelineGraph()
        BioToolkit.add_node!(pg, BioToolkit.pipeline_node(:n1, () -> 1))
        retried = BioToolkit.retry_execute_pipeline(pg; max_retries=1)
        @test retried.status == :ok

        slurm = BioToolkit.slurm_array_plan(DataFrame(sample_id=["s1", "s2", "s3"]); cpus_per_task=4, mem_gb=16)
        @test nrow(slurm) == 3
    end

    @testset "Somatic and perturbation" begin
        g = BioToolkit.germline_bayesian_call(18, 12)
        @test g.genotype in ("0/0", "0/1", "1/1")

        s = BioToolkit.somatic_tumor_normal_call(10, 15, 30, 0)
        @test s.somatic == true

        cn_s = BioToolkit.somatic_tumor_normal_call(75, 25, 100, 0; tumor_copy_number=4.0, tumor_purity=0.5)
        @test cn_s.tumor_expected_vaf_somatic ≈ 1 / 6 atol=1e-6
        @test cn_s.posterior_somatic >= 0.5

        hc = BioToolkit.haplotypecaller_call("ref.fa", "normal.bam")
        @test hc.status == :planned

        mt = BioToolkit.mutect2_call("ref.fa", "tumor.bam"; normal_bam="normal.bam")
        @test mt.status == :planned

        st = BioToolkit.strelka2_call("ref.fa", "tumor.bam", "normal.bam")
        @test st.status == :planned

        sv = BioToolkit.sv_breakpoint_graph_table(DataFrame(chrom1=["chr1", "chr1"], pos1=[100, 100], chrom2=["chr2", "chr2"], pos2=[200, 200]))
        @test nrow(sv) == 1

        ann = BioToolkit.vep_like_annotation(DataFrame(consequence=["missense_variant", "intron_variant"]))
        @test "impact_label" in names(ann)

        mix = BioToolkit.mixscape_like_contrast(rand(8, 10), ["control", "control", "drug", "drug", "drug", "control", "drug", "control", "drug", "control"])
        @test nrow(mix) > 0

        mixb = BioToolkit.mixscape_multibatch_contrast(rand(8, 12), ["control", "drug", "drug", "control", "drug", "control", "drug", "control", "drug", "control", "drug", "control"], ["b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2", "b3", "b3", "b3", "b3"])
        @test nrow(mixb) > 0

        pb = BioToolkit.perturbseq_pseudobulk(abs.(randn(6, 30)), ["g1", "g2", "g3", "g4", "g5", "g6"], vcat(fill("c", 10), fill("a", 10), fill("b", 10)); min_cells=5)
        @test size(pb.counts, 2) == 3

        sce2 = BioToolkit.SingleCellExperiment(Int.(round.(abs.(randn(12, 24)) .* 10)), ["g$(i)" for i in 1:12], ["c$(i)" for i in 1:24])
        milo = BioToolkit.milo_like_neighborhood_da(sce2, vcat(fill("ctrl", 12), fill("case", 12)); k=5, reference_group="ctrl")
        @test nrow(milo) == 24

        syn = BioToolkit.perturbation_synergy_scores(rand(10, 40), vcat(fill(false, 10), fill(true, 10), fill(false, 10), fill(true, 10)), vcat(fill(false, 10), fill(false, 10), fill(true, 10), fill(true, 10)))
        @test nrow(syn) == 10

        cg = BioToolkit.causal_grn_inference(rand(12, 20); top_k=30)
        @test nrow(cg) <= 30

        nmf = BioToolkit.gene_program_nmf(abs.(randn(25, 14)); n_programs=4, n_iter=50, backend=:cpu)
        @test size(nmf.program_loadings, 2) == 4

        tmb = BioToolkit.tumor_mutational_burden(DataFrame(sample=["a", "a", "b"], consequence=["missense_variant", "intron_variant", "frameshift_variant"]); callable_mb=30.0)
        @test nrow(tmb) == 2

        sig = BioToolkit.mutational_signature_nmf(abs.(randn(96, 5)); n_signatures=3, n_iter=40, backend=:cpu)
        @test size(sig.signatures, 2) == 3

        seg = BioToolkit.cbs_like_segments(randn(80); min_bins=8)
        @test nrow(seg) >= 1
    end
end
