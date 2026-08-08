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
end
