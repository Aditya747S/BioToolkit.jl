using Test
using SparseArrays
using DataFrames

@testset "Repertoire" begin
    @testset "RepertoireContig and ContigData types" begin
        contig = BioToolkit.RepertoireContig(
            "contig1", "cell1", :TRB, "TRBV12-3*01", "TRBD1*01", "TRBJ2-1*01", "TRBC1*01",
            "GATGCT", "GCCAGC", "ATGGCC",
            "DA", "AS", "CAVRDGADHTDTQYF",
            "ATGGCCGTGAGAGACGGTGCTGACCACACCCAGTACTTT", 40, 5, 10, 0.95, true, true,
            Dict{String,Any}()
        )
        @test contig.contig_id == "contig1"
        @test contig.barcode == "cell1"
        @test contig.chain == :TRB
        @test contig.v_gene == "TRBV12-3*01"
        @test contig.cdr3_aa == "CAVRDGADHTDTQYF"
        @test contig.length == 40
        @test contig.umi_count == 5
        @test contig.productive == true

        contig2 = BioToolkit.RepertoireContig(
            "contig2", "cell1", :TRA, "TRAV12-1*01", "", "TRAJ26*01", "TRAC1*01",
            "", "", "ATGGCC",
            "", "", "CAVSGYNQGGTSGCSYTLTF",
            "ATGGCCGTGAGCGGTTACAACGGCGGCACCTCGGGCTGCTCTTACACTTTGACCTTT", 52, 3, 8, 0.88, true, false,
            Dict{String,Any}()
        )

        contigs = [contig, contig2]
        data = BioToolkit.ContigData(contigs, "sample1")
        @test length(data) == 2
        @test !isempty(data)
        @test data.sample_id == "sample1"
    end

    @testset "AIRR parsing" begin
        airr_row = Dict(
            "sequence_id" => "seq1",
            "cell_id" => "cell1",
            "locus" => "TRA",
            "v_call" => "TRAV12-1*01",
            "d_call" => "",
            "j_call" => "TRAJ26*01",
            "c_call" => "TRAC1*01",
            "junction_aa" => "CAVSGYNQGGTSGCSYTLTF",
            "sequence" => "atggccgtgagcggttacaacggcggcacctcgggctgctcttacactttgaccttt",
            "duplicate_count" => 3,
            "productive" => true
        )
        parsed = BioToolkit.airr_parse_row(airr_row)
        @test parsed.barcode == "cell1"
        @test parsed.chain == :TRA
        @test parsed.cdr3_aa == "CAVSGYNQGGTSGCSYTLTF"
        @test parsed.umi_count == 3
        @test parsed.productive == true
    end

    @testset "Diversity metrics" begin
        clone_sizes = [100, 50, 30, 20, 10, 5, 3, 2, 1, 1]
        total = sum(clone_sizes)

        h = BioToolkit.shannon_entropy(clone_sizes)
        @test h > 0.0
        @test h <= log(length(clone_sizes))

        gs = BioToolkit.gini_simpson(clone_sizes)
        @test 0.0 <= gs <= 1.0

        clonality = BioToolkit.clonality_index(clone_sizes)
        @test 0.0 <= clonality <= 1.0

        chao1 = BioToolkit.chao1_richness(clone_sizes)
        @test chao1 >= length(clone_sizes)

        ace = BioToolkit.ace_richness(clone_sizes)
        @test ace >= length(clone_sizes)

        hill0 = BioToolkit.hill_diversity(clone_sizes, 0.0)
        hill1 = BioToolkit.hill_diversity(clone_sizes, 1.0)
        hill2 = BioToolkit.hill_diversity(clone_sizes, 2.0)
        @test hill0 == length(clone_sizes)
        @test hill1 > 0.0
        @test hill2 > 0.0
        @test hill1 <= exp(h)
        @test hill2 <= total^2 / minimum(clone_sizes)^2

        diversity_df = BioToolkit.diversity_curve(clone_sizes; q_values=[0.0, 1.0, 2.0, 3.0])
        @test nrow(diversity_df) == 4
        @test diversity_df.q[1] == 0.0
        @test diversity_df.q[2] == 1.0

        rare = BioToolkit.rarefaction_curve(clone_sizes; n_points=10)
        @test hasproperty(rare, :sample_size)
        @test hasproperty(rare, :rarefaction)
        @test hasproperty(rare, :extrapolation)
        @test all(rare.rarefaction .>= 0)
        @test all(rare.extrapolation .>= 0)

        singletons = [1, 1, 1, 1, 1]
        h_single = BioToolkit.shannon_entropy(singletons)
        gs_single = BioToolkit.gini_simpson(singletons)
        @test h_single ≈ log(5)
        @test gs_single ≈ 0.8

        uniform = [10, 10, 10, 10]
        h_uniform = BioToolkit.shannon_entropy(uniform)
        gs_uniform = BioToolkit.gini_simpson(uniform)
        @test h_uniform ≈ log(4)
        @test gs_uniform ≈ 0.75
    end

    @testset "Clonotype combination" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "TRBV12-3*01", "", "TRBJ2-1*01", "", "", "", "", "", "", "CVRGGGADHTQYF", "ATGTGTCGAGGGGGCGCTGACCACACCCAGTAC", 34, 5, 10, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell1", :TRA, "TRAV12-1*01", "", "TRAJ26*01", "", "", "", "", "", "", "CAVSGYNQGGTSGCSYTLTF", "ATG", 3, 2, 5, 0.8, true, false, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell2", :TRB, "TRBV12-3*01", "", "TRBJ2-1*01", "", "", "", "", "", "", "CVRTYYGSSYEQYF", "ATGTGTCGTACTTATGGGAGCTCCTATGAGCAATTT", 38, 8, 15, 0.95, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell2", :TRA, "TRAV12-1*01", "", "TRAJ26*01", "", "", "", "", "", "", "CAVSGYNQGGTSGCSYTLTF", "ATG", 3, 1, 3, 0.7, true, false, Dict{String,Any}()),
        ]
        data = BioToolkit.ContigData(contigs, "sample1")
        combined = BioToolkit.combine_clonotypes(data; chain_pair=:TRB)
        @test length(combined) == 2
    end

    @testset "CDR3 Levenshtein distance" begin
        seqs = ["CAVSGYNQGGTSGCSYTLTF", "CAVSGANQGGTSGCSYTLTF", "CAVSGYNQGGTSG", "TRAVSGYNQGGTSGCSY"]
        D = BioToolkit.cdr3_levenshtein_distance(seqs)
        @test size(D) == (4, 4)
        @test D[1, 1] == 0.0
        @test D[1, 2] == 1.0
        @test D[1, 3] > D[1, 2]
        @test D == D'
    end

    @testset "CDR3 clustering" begin
        seqs = ["CVVGGG", "CVGGG", "CAAGGG", "TRAAGGG", "CVVGGG", "CVGGG", "TRBVGGG"]
        result = BioToolkit.cdr3_clustering(seqs; threshold=2, min_size=2)
        @test haskey(result, :assignments)
        @test haskey(result, :clusters)
        @test length(result.assignments) == 7
    end

    @testset "Clone size distribution" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 8, 15, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 5, 10, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATGTT", 5, 3, 6, 0.85, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c5", "cell5", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATGTT", 5, 2, 4, 0.85, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c6", "cell6", :TRB, "V3", "", "J3", "", "", "", "", "", "", "C3", "ATGCAG", 6, 1, 2, 0.8, true, true, Dict{String,Any}()),
        ]
        data = BioToolkit.ContigData(contigs, "sample1")
        sample = BioToolkit.RepertoireSample(data)

        @test length(sample) == 6
        @test !isempty(sample.clonotypes)
        @test hasproperty(sample.clonotypes, :clone_size)

        summary = BioToolkit.clonal_distribution_summary(sample)
        @test hasproperty(summary, :metric)
        @test hasproperty(summary, :value)
        @test nrow(summary) > 0

        stats = BioToolkit.repertoire_statistics(sample)
        @test hasproperty(stats, :shannon_entropy)
        @test hasproperty(stats, :gini_simpson)
        @test hasproperty(stats, :clonality)
        @test stats[1, :n_cells] == 6
        @test stats[1, :n_clonotypes] == 3

        expanded = BioToolkit.expanded_clonotype_test(sample; threshold=5)
        @test hasproperty(expanded, :pvalue)
        @test hasproperty(expanded, :padj)
    end

    @testset "V(J) usage analysis" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "TRBV12-3*01", "", "TRBJ2-1*01", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "TRBV12-3*01", "", "TRBJ2-1*01", "", "", "", "", "", "", "C1", "ATGGCC", 6, 8, 15, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "TRBV12-5*01", "", "TRBJ2-3*01", "", "", "", "", "", "", "C2", "ATGTT", 5, 5, 10, 0.85, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "TRBV12-5*01", "", "TRBJ2-3*01", "", "", "", "", "", "", "C2", "ATGTT", 5, 3, 6, 0.85, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c5", "cell5", :TRB, "TRBV12-5*01", "", "TRBJ2-1*01", "", "", "", "", "", "", "C3", "ATGCAG", 6, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]
        data = BioToolkit.ContigData(contigs, "sample1")
        sample = BioToolkit.RepertoireSample(data)

        vj = BioToolkit.vj_usage_matrix(sample)
        @test haskey(vj, :v_usage)
        @test haskey(vj, :j_usage)
        @test nrow(vj.v_usage) > 0
        @test hasproperty(vj.v_usage, :frequency)

        pairs = BioToolkit.vj_pairing_analysis(sample)
        @test hasproperty(pairs, :v_gene)
        @test hasproperty(pairs, :j_gene)
        @test hasproperty(pairs, :count)
    end

    @testset "Overlap metrics" begin
        contigs1 = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 8, 15, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATGTT", 5, 5, 10, 0.85, true, true, Dict{String,Any}()),
        ]
        contigs2 = [
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c5", "cell5", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATGTT", 5, 5, 10, 0.85, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c6", "cell6", :TRB, "V3", "", "J3", "", "", "", "", "", "", "C3", "ATGCAG", 6, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]

        data1 = BioToolkit.ContigData(contigs1, "sample1")
        data2 = BioToolkit.ContigData(contigs2, "sample2")
        sample1 = BioToolkit.RepertoireSample(data1)
        sample2 = BioToolkit.RepertoireSample(data2)

        set_a = Set{String}(sample1.clonotypes.clonotype_id)
        set_b = Set{String}(sample2.clonotypes.clonotype_id)
        jaccard = BioToolkit.overlap_jaccard(set_a, set_b)
        @test 0.0 <= jaccard <= 1.0
        @test jaccard > 0.0

        mh = BioToolkit.overlap_morisita_horn(sample1, sample2)
        @test 0.0 <= mh <= 1.0

        dist_jaccard = BioToolkit.repertoire_distance(sample1, sample2; metric=:jaccard)
        @test 0.0 <= dist_jaccard <= 1.0

        samples = [sample1, sample2]
        overlap_mat = BioToolkit.clonotype_overlap_matrix(samples; metric=:jaccard)
        @test size(overlap_mat) == (2, 2)
        @test overlap_mat[1, 1] == 1.0
        @test overlap_mat[1, 2] == overlap_mat[2, 1]

        divergence = BioToolkit.repertoire_divergence(samples; metric=:jaccard)
        @test size(divergence) == (2, 2)
    end

    @testset "Convergent clonotypes" begin
        contigs1 = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CV1J1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CV1J1", "ATGGCC", 6, 8, 15, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V2", "", "J2", "", "", "", "", "", "", "CV2J2", "ATGTT", 5, 5, 10, 0.85, true, true, Dict{String,Any}()),
        ]
        contigs2 = [
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CV1J1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c5", "cell5", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CV1J1", "ATGGCC", 6, 7, 12, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c6", "cell6", :TRB, "V3", "", "J3", "", "", "", "", "", "", "CV3J3", "ATGCAG", 6, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]

        data1 = BioToolkit.ContigData(contigs1, "sample1")
        data2 = BioToolkit.ContigData(contigs2, "sample2")
        sample1 = BioToolkit.RepertoireSample(data1)
        sample2 = BioToolkit.RepertoireSample(data2)

        convergent = BioToolkit.convergent_clonotypes([sample1, sample2]; min_shared=1)
        @test hasproperty(convergent, :clonotype_id)
        @test hasproperty(convergent, :n_samples)
        @test any(n >= 2 for n in convergent.n_samples)
    end

    @testset "Advanced Overlap & Divergence Metrics" begin
        contigs1 = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATGTT", 5, 5, 10, 0.85, true, true, Dict{String,Any}()),
        ]
        contigs2 = [
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATGGCC", 6, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "V3", "", "J3", "", "", "", "", "", "", "C3", "ATGCAG", 6, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]
        sample1 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs1, "sample1"))
        sample2 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs2, "sample2"))

        set_a = Set{String}(sample1.clonotypes.clonotype_id)
        set_b = Set{String}(sample2.clonotypes.clonotype_id)

        dice = BioToolkit.overlap_sorensen_dice(set_a, set_b)
        @test 0.0 <= dice <= 1.0

        oc = BioToolkit.overlap_overlap_coefficient(set_a, set_b)
        @test 0.0 <= oc <= 1.0

        cos_sim = BioToolkit.overlap_cosine(sample1, sample2)
        @test 0.0 <= cos_sim <= 1.0

        heatmap_df = BioToolkit.repertoire_heatmap([sample1, sample2])
        @test nrow(heatmap_df) == 4
        @test hasproperty(heatmap_df, :value)

        venn_df = BioToolkit.clonotype_venn([sample1, sample2])
        @test nrow(venn_df) > 0
        @test hasproperty(venn_df, :count)
    end

    @testset "TCRdist & Hamming Matrix" begin
        seqs = ["CAVSGYNQGGTSGCSYTLTF", "CAVSGANQGGTSGCSYTLTF", "CAVSGYNQGGTSGCSYTLTT"]
        H = BioToolkit.cdr3_hamming_distance(seqs)
        @test size(H) == (3, 3)
        @test H[1, 2] == 1.0

        T = BioToolkit.tcrdist_matrix(seqs)
        @test size(T) == (3, 3)
        @test T[1, 1] == 0.0
        @test T[1, 2] > 0.0
    end

    @testset "Clonotype Network Graph" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CAVSGYNQGGTSGCSYTLTF", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CAVSGANQGGTSGCSYTLTF", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V2", "", "J2", "", "", "", "", "", "", "CVVGGGGGGGGGGGGGGGGG", "ATG", 3, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]
        sample = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs, "sample1"))
        net = BioToolkit.clonotype_network(sample; threshold=2)

        @test haskey(net, :adjacency)
        @test haskey(net, :nodes)
        @test haskey(net, :edges)
        @test haskey(net, :components)
        @test nrow(net.nodes) == 3
    end

    @testset "Public Clonotype & VJ Divergence" begin
        contigs1 = [BioToolkit.RepertoireContig("c1", "cell1", :TRB, "TRBV12-1*01", "", "TRBJ1-1*01", "", "", "", "", "", "", "C1", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}())]
        contigs2 = [BioToolkit.RepertoireContig("c2", "cell2", :TRB, "TRBV12-1*01", "", "TRBJ1-1*01", "", "", "", "", "", "", "C1", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}())]
        sample1 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs1, "sample1"))
        sample2 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs2, "sample2"))

        pub = BioToolkit.public_clonotype_analysis([sample1, sample2]; min_samples=2)
        @test nrow(pub) == 1
        @test pub.prevalence[1] == 2

        jsd = BioToolkit.vj_usage_divergence(sample1, sample2; metric=:jsd)
        @test jsd >= 0.0

        diff_vj = BioToolkit.differential_vj_usage([sample1], [sample2])
        @test hasproperty(diff_vj, :v_gene)
        @test hasproperty(diff_vj, :pvalue)
    end

    @testset "BCR SHM Rate & Isotype Summary" begin
        bcr_contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :IGH, "IGHV1-2*01", "IGHD2-2*01", "IGHJ4*01", "IGHM", "", "", "", "", "", "CAR", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}("shm_rate" => 0.03)),
            BioToolkit.RepertoireContig("c2", "cell2", :IGH, "IGHV1-2*01", "IGHD2-2*01", "IGHJ4*01", "IGHG1", "", "", "", "", "", "CAR", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}("shm_rate" => 0.05)),
        ]
        shm_res = BioToolkit.bcr_shm_rate(bcr_contigs)
        @test shm_res.mean_shm > 0.0
        @test nrow(shm_res.shm_df) == 2

        iso_res = BioToolkit.isotype_usage_summary(bcr_contigs)
        @test nrow(iso_res) == 2
        @test hasproperty(iso_res, :isotype)
        @test hasproperty(iso_res, :frequency)
    end

    @testset "BKTree Metric Search" begin
        seqs = ["CAVSGYNQGGTSGCSYTLTF", "CAVSGANQGGTSGCSYTLTF", "CAVSGYNQGGTSG", "CVVGGGGGGGGGGGGGGGGG"]
        tree = BioToolkit.BKTree(seqs)
        res = BioToolkit.bktree_range_query(tree, "CAVSGYNQGGTSGCSYTLTF", 2)
        @test 1 in res
        @test 2 in res
        @test !(4 in res)
    end

    @testset "Multi-Chain Pairing Analysis" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRA, "TRAV1", "", "TRAJ1", "", "", "", "", "", "", "C1", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell1", :TRB, "TRBV1", "", "TRBJ1", "", "", "", "", "", "", "C2", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c3", "cell2", :TRA, "TRAV2", "", "TRAJ2", "", "", "", "", "", "", "C3", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell2", :TRA, "TRAV3", "", "TRAJ3", "", "", "", "", "", "", "C4", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c5", "cell2", :TRB, "TRBV2", "", "TRBJ2", "", "", "", "", "", "", "C5", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
        ]
        data = BioToolkit.ContigData(contigs, "sample1")
        pairing = BioToolkit.multi_chain_pairing_analysis(data)
        @test nrow(pairing) == 5
        @test pairing.cell_count[1] == 1 # Canonical
        @test pairing.cell_count[2] == 1 # Dual Alpha
    end

    @testset "Bootstrap Diversity CIs & Permutation Overlap Test" begin
        contigs1 = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V2", "", "J2", "", "", "", "", "", "", "C2", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
        ]
        contigs2 = [
            BioToolkit.RepertoireContig("c3", "cell3", :TRB, "V1", "", "J1", "", "", "", "", "", "", "C1", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c4", "cell4", :TRB, "V3", "", "J3", "", "", "", "", "", "", "C3", "ATG", 3, 2, 4, 0.8, true, true, Dict{String,Any}()),
        ]
        sample1 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs1, "sample1"))
        sample2 = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs2, "sample2"))

        ci_df = BioToolkit.bootstrap_diversity_ci(sample1; n_bootstraps=50)
        @test nrow(ci_df) == 4
        @test hasproperty(ci_df, :ci_lower)
        @test hasproperty(ci_df, :ci_upper)

        perm_test = BioToolkit.repertoire_overlap_permutation_test(sample1, sample2; n_permutations=100)
        @test haskey(perm_test, :observed_overlap)
        @test haskey(perm_test, :pvalue)
        @test 0.0 <= perm_test.pvalue <= 1.0
    end

    @testset "GLIPH2 Motif Enrichment & Positional Logos" begin
        cdr3s = ["CAVSGYNQGGTSGCSYTLTF", "CAVSGANQGGTSGCSYTLTF", "CAVSGYNQGGTSGCSYTLTT", "CAVSGYNQGGTSGCSYTLTA"]
        motifs = BioToolkit.gliph2_motif_enrichment(cdr3s; k_range=3:3, p_cutoff=1.0)
        @test nrow(motifs) > 0
        @test hasproperty(motifs, :motif)
        @test hasproperty(motifs, :pvalue)

        logo = BioToolkit.positional_motif_logo(cdr3s)
        @test nrow(logo) == 20
        @test hasproperty(logo, :information_bits)
    end

    @testset "CDR3 Biophysical & Atchley Factors" begin
        cdr3s = ["CAVSGYNQGGTSGCSYTLTF", "CVVGGGGGGGGGGGGGGGGG"]
        props = BioToolkit.repertoire_cdr3_physicochemical_properties(cdr3s)
        @test nrow(props) == 2
        @test hasproperty(props, :gravy_hydrophobicity)
        @test hasproperty(props, :net_charge)
        @test hasproperty(props, :atchley_polarity)
        @test hasproperty(props, :atchley_charge)
    end

    @testset "Spectratyping & Lineage Trees" begin
        contigs = [
            BioToolkit.RepertoireContig("c1", "cell1", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CAVSGYNQGGTSGCSYTLTF", "ATG", 3, 10, 20, 0.9, true, true, Dict{String,Any}()),
            BioToolkit.RepertoireContig("c2", "cell2", :TRB, "V1", "", "J1", "", "", "", "", "", "", "CAVSGANQGGTSGCSYTLTF", "ATG", 3, 5, 10, 0.9, true, true, Dict{String,Any}()),
        ]
        sample = BioToolkit.RepertoireSample(BioToolkit.ContigData(contigs, "sample1"))

        spec = BioToolkit.spectratype_analysis(sample)
        @test haskey(spec, :spectratype)
        @test haskey(spec, :spectratype_entropy)
        @test spec.spectratype_entropy >= 0.0

        tree = BioToolkit.clonotype_lineage_tree(sample, sample.clonotypes.clonotype_id[1])
        @test haskey(tree, :tree_nodes)
        @test haskey(tree, :tree_edges)
    end

    @testset "BCR Mutation Hotspot Profile" begin
        seqs = ["AGCTACAAAGCTACAA", "CCCCCCCCCCCCCCCC"]
        hotspots = BioToolkit.bcr_mutation_hotspot_profile(seqs)
        @test nrow(hotspots) == 2
        @test hasproperty(hotspots, :hotspot_count)
        @test hotspots.hotspot_count[1] > hotspots.hotspot_count[2]
    end
end

