using Test
using SparseArrays
using DataFrames
using LinearAlgebra

@testset "Hi-C" begin
    @testset "HiCContactMatrix types" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", 0, 9999, '.'),
            BioToolkit.GenomicInterval("chr1", 10000, 19999, '.'),
            BioToolkit.GenomicInterval("chr1", 20000, 29999, '.'),
            BioToolkit.GenomicInterval("chr1", 30000, 39999, '.'),
        ]
        m = sparse([0 10 5 2; 10 0 8 3; 5 8 0 6; 2 3 6 0])
        matrix = BioToolkit.HiCContactMatrix(m, "chr1", bins, 10000, false, "none", Dict{String,Any}())
        @test length(matrix) == 4
        @test size(matrix) == (4, 4)
        @test matrix.chrom == "chr1"
        @test matrix.bin_size == 10000
        @test matrix.normalized == false
        @test length(matrix.bins) == 4
    end

    @testset "Normalization" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", 0, 9999, '.'),
            BioToolkit.GenomicInterval("chr1", 10000, 19999, '.'),
            BioToolkit.GenomicInterval("chr1", 20000, 29999, '.'),
            BioToolkit.GenomicInterval("chr1", 30000, 39999, '.'),
        ]
        m = sparse([0 10 5 2; 10 0 8 3; 5 8 0 6; 2 3 6 0])
        matrix = BioToolkit.HiCContactMatrix(m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        ice_norm = BioToolkit.normalize_ice(matrix)
        @test ice_norm isa BioToolkit.HiCContactMatrix
        @test ice_norm.normalized == true
        @test ice_norm.normalization_method == "ICE"
        @test size(ice_norm.matrix) == size(m)

        kr_norm = BioToolkit.normalize_kr(matrix)
        @test kr_norm isa BioToolkit.HiCContactMatrix
        @test kr_norm.normalized == true
        @test kr_norm.normalization_method == "KR"

        vanilla_norm = BioToolkit.normalize_vanilla(matrix)
        @test vanilla_norm isa BioToolkit.HiCContactMatrix
        @test vanilla_norm.normalized == true
        @test vanilla_norm.normalization_method == "vanilla"
    end

    @testset "Expected decay and power law" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:49
        ]
        n = length(bins)

        decay_m = spzeros(n, n)
        for i in 1:n
            for j in (i+1):n
                d = j - i
                decay_m[i, j] = 100.0 / (d^0.8)
                decay_m[j, i] = decay_m[i, j]
            end
        end
        decay_m[diagind(decay_m)] .= 1000.0

        matrix = BioToolkit.HiCContactMatrix(decay_m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        expected_df = BioToolkit.expected_decay(matrix)
        @test expected_df isa DataFrame
        @test hasproperty(expected_df, :distance_bins)
        @test hasproperty(expected_df, :genomic_distance)
        @test hasproperty(expected_df, :expected)
        @test nrow(expected_df) > 0
        @test all(expected_df.expected .>= 0)

        fit = BioToolkit.power_law_fit(expected_df)
        @test haskey(fit, :alpha)
        @test haskey(fit, :r_squared)
        @test isfinite(fit.alpha)
        @test isfinite(fit.r_squared)
        @test 0.0 <= fit.alpha <= 3.0

        oe = BioToolkit.obs_exp_ratio(matrix, expected_df)
        @test oe isa BioToolkit.HiCContactMatrix
        @test size(oe.matrix) == (n, n)
    end

    @testset "Insulation score and boundaries" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:19
        ]
        n = length(bins)

        contact_m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                contact_m[i, j] = max(1.0, 100.0 / (d + 1)^0.8)
            end
        end
        contact_m[diagind(contact_m)] .= 500.0

        matrix = BioToolkit.HiCContactMatrix(contact_m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        ins = BioToolkit.hic_insulation_score(matrix; window_bins=5)
        @test ins isa DataFrame
        @test hasproperty(ins, :bin)
        @test hasproperty(ins, :insulation)
        @test hasproperty(ins, :zscore)
        @test nrow(ins) == n

        ins = BioToolkit.hic_insulation_score(matrix; window_bins=3)
        @test all(isfinite, ins.insulation[isfinite.(ins.insulation)])

        boundaries = BioToolkit.boundary_detection(ins; q_threshold=0.2, min_strength=0.1)
        @test boundaries isa DataFrame
        @test hasproperty(boundaries, :position)
        @test hasproperty(boundaries, :strength)
    end

    @testset "Compartments" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:19
        ]
        n = length(bins)

        compartment_m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                if (i <= 10 && j <= 10) || (i > 10 && j > 10)
                    compartment_m[i, j] = 100.0 / (abs(i - j) + 1)^0.8
                else
                    compartment_m[i, j] = 5.0 / (abs(i - j) + 1)^0.8
                end
            end
        end

        matrix = BioToolkit.HiCContactMatrix(compartment_m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        scores = BioToolkit.compartment_score(matrix)
        @test scores isa DataFrame
        @test hasproperty(scores, :pc1)
        @test hasproperty(scores, :compartment)
        @test length(scores.pc1) == n
        @test all(c -> c in ("A", "B"), scores.compartment)

        compartments = BioToolkit.ab_compartments(matrix)
        @test haskey(compartments, :scores)
        @test haskey(compartments, :a_fraction)
        @test 0.0 <= compartments.a_fraction <= 1.0
        @test 0.0 <= compartments.b_fraction <= 1.0
        @test abs(compartments.a_fraction + compartments.b_fraction - 1.0) < 0.01

        compartments2 = BioToolkit.ab_compartments(matrix)
        switching = BioToolkit.compartment_switching(scores, compartments2.scores; min_flip=0.1)
        @test switching isa DataFrame
        @test hasproperty(switching, :bin)
        @test hasproperty(switching, :switch_type)
    end

    @testset "Loop calling" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:39
        ]
        n = length(bins)

        loop_m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                loop_m[i, j] = 100.0 / (d + 1)^0.8
            end
        end
        loop_m[10, 30] = 1000.0
        loop_m[30, 10] = 1000.0
        loop_m[15, 35] = 800.0
        loop_m[35, 15] = 800.0

        matrix = BioToolkit.HiCContactMatrix(loop_m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        loops = BioToolkit.call_loops(matrix; window=5, fdr=0.5, min_dist=3, max_dist=20)
        @test loops isa DataFrame
        @test hasproperty(loops, :anchor1_bin)
        @test hasproperty(loops, :anchor2_bin)
        @test hasproperty(loops, :score)
        @test all(loops.score .>= 0)

        loops2 = BioToolkit.call_loops(matrix; window=5, fdr=0.5, min_dist=3, max_dist=20)
        merged = BioToolkit.merge_loops(loops, loops2; max_distance=3)
        @test merged isa DataFrame
        @test hasproperty(merged, :found_in)
        @test hasproperty(merged, :max_score)
    end

    @testset "Pileup and APA" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:29
        ]
        n = length(bins)

        pileup_m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                pileup_m[i, j] = 100.0 / (d + 1)^0.8
            end
        end
        pileup_m[15, 15] = 500.0

        matrix = BioToolkit.HiCContactMatrix(pileup_m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        loops_df = DataFrame(
            anchor1_bin=[15, 20],
            anchor2_bin=[15, 20],
            anchor1_pos=[150000, 200000],
            anchor2_pos=[150000, 200000],
            score=[10.0, 8.0]
        )

        pileup = BioToolkit.pileup_matrix(matrix, loops_df; window_bins=5)
        @test haskey(pileup, :matrix)
        @test haskey(pileup, :n_loops)
        @test pileup.n_loops >= 0

        aparesult = BioToolkit.loop_apa(pileup)
        @test haskey(aparesult, :apa_score)
        @test haskey(aparesult, :center_mean)
        @test haskey(aparesult, :corner_mean)
        @test aparesult.apa_score >= 0.0
    end

    @testset "Differential interaction testing" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:19
        ]
        n = length(bins)

        m_a = spzeros(n, n)
        m_b = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                m_a[i, j] = max(1.0, 50.0 / (d + 1)^0.8)
                m_b[i, j] = max(1.0, 80.0 / (d + 1)^0.8)
            end
        end
        m_a[10, 10] = 200.0
        m_b[10, 10] = 1000.0
        m_a[10, 10] = 200.0
        m_b[10, 10] = 1000.0

        matrix_a = BioToolkit.HiCContactMatrix(m_a, "chr1", bins, 10000, false, "none", Dict{String,Any}())
        matrix_b = BioToolkit.HiCContactMatrix(m_b, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        coords = [(10, 10), (5, 5), (15, 15)]
        diff = BioToolkit.differential_interaction([matrix_a], [matrix_b], coords)
        @test diff isa DataFrame
        @test hasproperty(diff, :bin1)
        @test hasproperty(diff, :bin2)
        @test hasproperty(diff, :log2_fc)
        @test hasproperty(diff, :pvalue)
        @test hasproperty(diff, :padj)
        @test nrow(diff) == 3
    end

    @testset "Quality metrics" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:19
        ]
        n = length(bins)

        m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                m[i, j] = max(1.0, 100.0 / (d + 1)^0.8)
            end
        end

        matrix = BioToolkit.HiCContactMatrix(m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        qc = BioToolkit.hic_quality_metrics(matrix)
        @test qc isa DataFrame
        @test hasproperty(qc, :metric)
        @test hasproperty(qc, :value)
        @test any(qc.metric .== "total_contacts")
        @test any(qc.metric .== "n_bins")
        @test any(qc.metric .== "sparsity")

        qc = BioToolkit.hic_quality_metrics(matrix)
        @test qc isa DataFrame
    end

    @testset "Contact decay" begin
        bins = [
            BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:29
        ]
        n = length(bins)

        m = spzeros(n, n)
        for i in 1:n
            for j in 1:n
                d = abs(j - i)
                m[i, j] = max(1.0, 100.0 / (d + 1)^0.9)
            end
        end

        matrix = BioToolkit.HiCContactMatrix(m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

        decay = BioToolkit.contact_decay_by_distance(matrix)
        @test decay isa DataFrame
        @test hasproperty(decay, :expected)

        decay_limited = BioToolkit.contact_decay_by_distance(matrix; max_dist=10)
        @test nrow(decay_limited) <= nrow(decay)
    end

    @testset "Cool file operations" begin
        mktempdir() do dir
            bins = [
                BioToolkit.GenomicInterval("chr1", i * 10000, (i + 1) * 10000 - 1, '.') for i in 0:9
            ]
            n = length(bins)
            m = sparse([0 10 5 0 0 0 0 0 0 0;
                        10 0 8 0 0 0 0 0 0 0;
                        5 8 0 6 0 0 0 0 0 0;
                        0 0 6 0 4 0 0 0 0 0;
                        0 0 0 4 0 3 0 0 0 0;
                        0 0 0 0 3 0 2 0 0 0;
                        0 0 0 0 0 2 0 1 0 0;
                        0 0 0 0 0 0 1 0 5 0;
                        0 0 0 0 0 0 0 5 0 8;
                        0 0 0 0 0 0 0 0 8 0])
            matrix = BioToolkit.HiCContactMatrix(m, "chr1", bins, 10000, false, "none", Dict{String,Any}())

            cool_path = joinpath(dir, "test.cool")
            BioToolkit.write_cool(matrix, cool_path)
            @test isfile(cool_path)

            content = read(cool_path, String)
            @test occursin("chrom", content)
            @test occursin("count", content)
        end
    end
end
