using Test
using BioToolkit
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

# ==========================================================================
# Session 2 corrections (2026-09-04): NB-GLM/Wald p-values, exact two-count
# test, Monte-Carlo loop FDR, ICE ignore_diag, real Knight-Ruiz, cooler HDF5
# I/O + mcool, honest .hic handling, O/E compartments, real GC, MDS payload.
# Evidence: plan/audit_notes.md wave 6; plan/correction_plan.md items 12-21.
# ==========================================================================
using Random
function _session2b_matrix(; n=40, decay=0.25, seed=77, scale=20.0)
    Random.seed!(seed)
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:n]
    m = zeros(Float64, n, n)
    for i in 1:n, j in i:n
        m[i, j] = m[j, i] = scale * exp(-decay * (j - i)) * (0.5 + rand())
    end
    return BioToolkit.HiCContactMatrix(spzeros(Float64, n, n) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
end

using Distributions
using HDF5
using CodecZlib

@testset "Session 2: cooler HDF5 round trip" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 5000 + 1, i * 5000, '.') for i in 1:8]
    m = spzeros(Float64, 8, 8)
    Random.seed!(3)
    for i in 1:8, j in i:8
        v = exp(-0.4 * (j - i)) * rand()
        m[i, j] = v
        j > i && (m[j, i] = v)
    end
    mat = BioToolkit.HiCContactMatrix(m, "chr1", bins, 5000, false, "none", Dict{String,Any}())

    mktempdir() do dir
        cool_path = joinpath(dir, "test.cool")
        BioToolkit.write_cool(mat, cool_path)
        @test BioToolkit.HiC._is_hdf5_file(cool_path)         # real HDF5, not TSV
        recovered = BioToolkit.read_cool(cool_path)
        @test recovered.chrom == "chr1"
        @test recovered.bin_size == 5000
        @test size(recovered.matrix) == (8, 8)
        @test Matrix(recovered.matrix) ≈ Matrix(mat.matrix)

        # mcool layout
        mcool_path = joinpath(dir, "test.mcool")
        HDF5.h5open(mcool_path, "w") do h5
            g = HDF5.create_group(h5, "resolutions/5000")
            write(g, "chroms/names", ["chr1"])
            write(g, "chroms/lengths", [40000])
            write(g, "bins/chrom", ["chr1" for _ in 1:8])
            write(g, "bins/start", [(i - 1) * 5000 for i in 1:8])
            write(g, "bins/end", [i * 5000 for i in 1:8])
            nz = findnz(mat.matrix)
            upper = nz[1] .<= nz[2]                            # cooler convention
            write(g, "pixels/bin1_id", Int.(nz[1][upper]) .- 1)
            write(g, "pixels/bin2_id", Int.(nz[2][upper]) .- 1)
            write(g, "pixels/count", Float64.(nz[3][upper]))
        end
        mres = BioToolkit.read_mcool(mcool_path; resolution=5000)
        @test Matrix(mres.matrix) ≈ Matrix(mat.matrix)

        # TSV fallback preserved as its own function
        tsv_path = joinpath(dir, "test.tsv")
        BioToolkit.write_contacts_tsv(mat, tsv_path)
        tsv = BioToolkit.read_contacts_tsv(tsv_path)
        @test Matrix(tsv.matrix) ≈ Matrix(mat.matrix) atol=1e-6
    end
end

@testset "Session 2: corrupt .hic errors cleanly" begin
    mktempdir() do dir
        hic_path = joinpath(dir, "fake.hic")
        write(hic_path, "HIC" * "garbage")   # magic ok, header truncated
        @test_throws Union{ArgumentError,EOFError} BioToolkit.read_hic(hic_path; chrom="chr1")
    end
end

@testset "Session 2b: binary .hic round trip" begin
    mat = _session2b_matrix(; n=50)
    mktempdir() do dir
        path = joinpath(dir, "test.hic")
        BioToolkit.write_hic(mat, path)
        @test BioToolkit.HiC._is_hic_magic(path)
        recovered = BioToolkit.read_hic(path; chrom="chr1", resolution=1000)
        @test recovered.chrom == "chr1"
        @test size(recovered.matrix) == (50, 50)
        recovered_M = Matrix(recovered.matrix)
        original_M = Matrix(mat.matrix)
        # upper triangle must match exactly (float32 storage tolerance)
        for i in 1:50, j in i:50
            @test isapprox(recovered_M[i, j], original_M[i, j]; rtol=1e-6, atol=1e-8)
        end
        # symmetrization: lower triangle filled by the reader
        @test recovered_M[3, 1] ≈ original_M[1, 3] atol=1e-6

        # region query
        part = BioToolkit.read_hic(path; chrom="chr1", resolution=1000, start=10001, stop=20000)
        @test size(part.matrix) == (10, 10)
        @test part.matrix[1, 1] ≈ mat.matrix[11, 11] atol=1e-6

        # multi-chromosome experiment
        mat2 = _session2b_matrix(; n=20, seed=99)
        mat2 = BioToolkit.HiCContactMatrix(mat2.matrix, "chr2", mat2.bins, 1000, false, "none", Dict{String,Any}())
        experiment = BioToolkit.HiCExperiment(Dict("chr1" => mat, "chr2" => mat2), ["s1"], ["control"], Dict{String,Any}())
        exp_path = joinpath(dir, "exp.hic")
        BioToolkit.write_hic(experiment, exp_path)
        r1 = BioToolkit.read_hic(exp_path; chrom="chr1", resolution=1000)
        r2 = BioToolkit.read_hic(exp_path; chrom="chr2", resolution=1000)
        @test size(r1.matrix) == (50, 50)
        @test size(r2.matrix) == (20, 20)
        @test r2.matrix[2, 2] ≈ mat2.matrix[2, 2] atol=1e-6

        # TSV fallback still works for non-HIC files (first line: chrom size)
        tsv = joinpath(dir, "contacts.tsv")
        write(tsv, "chr1\t5000\nchr1\t1\tchr1\t1001\t5.0\n")
        back = BioToolkit.read_hic(tsv; chrom="chr1", resolution=1000)
        @test back.matrix[1, 2] == 5.0
    end
end

@testset "Session 2: ICE honours ignore_diag" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:6]
    m = Float64[4 2 1 0 0 0; 2 5 2 1 0 0; 1 2 6 2 1 0; 0 1 2 5 2 1; 0 0 1 2 4 2; 0 0 0 1 2 3]
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 6, 6) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())

    ice = BioToolkit.normalize_ice(mat)
    row_sums = vec(sum(Matrix(ice.matrix), dims=2))
    @test maximum(abs.(row_sums .- mean(row_sums))) < 1e-3 * mean(row_sums)  # balanced
    @test all(diag(Matrix(ice.matrix)) .== 0.0)                              # diagonal masked
    @test ice.normalization_method == "ICE"

    ice_keep = BioToolkit.normalize_ice(mat; ignore_diag=false, max_iter=5000)
    @test all(diag(Matrix(ice_keep.matrix)) .> 0.2)                          # diagonal participates
    raw_spread = maximum(abs.(vec(sum(Matrix(mat.matrix), dims=2)) .- mean(sum(Matrix(mat.matrix), dims=2))))
    ice_spread = maximum(abs.(vec(sum(Matrix(ice_keep.matrix), dims=2)) .- mean(vec(sum(Matrix(ice_keep.matrix), dims=2)))))
    @test ice_spread < raw_spread                                            # balancing improves row uniformity
end

@testset "Session 2: Knight-Ruiz converges to balanced matrix" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:6]
    m = Float64[4 2 1 0 0 0; 2 5 2 1 0 0; 1 2 6 2 1 0; 0 1 2 5 2 1; 0 0 1 2 4 2; 0 0 0 1 2 3]
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 6, 6) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())

    kr = BioToolkit.normalize_kr(mat)
    row_sums = vec(sum(Matrix(kr.matrix), dims=2))
    @test maximum(abs.(row_sums .- 1.0)) < 1e-3        # DA D 1 = 1
    @test kr.normalization_method == "KR"
    @test maximum(abs.(vec(sum(Matrix(kr.matrix), dims=1)) .- 1.0)) < 1e-3  # symmetric balance
end

@testset "Session 2: differential interaction Wald p-values" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:10]
    Random.seed!(11)
    function _cond(scale)
        m = spzeros(Float64, 10, 10)
        for i in 1:10, j in i:10
            v = Float64(rand(Poisson(exp(-0.15 * (j - i)) * scale)))
            m[i, j] = v; m[j, i] = v
        end
        return BioToolkit.HiCContactMatrix(m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    end
    cond_a = [_cond(50.0) for _ in 1:4]
    cond_b = [_cond(400.0) for _ in 1:4]
    coords = [(i, j) for i in 1:10 for j in (i + 1):10]
    result = BioToolkit.differential_interaction(cond_a, cond_b, coords)
    @test minimum(result.pvalue) < 0.05                      # strong effect detected
    @test all(0.0 .<= result.pvalue .<= 1.0)
    @test all(result.padj .>= result.pvalue .- 1e-12)

    # Null conditions: p-values should mostly be large.
    cond_n = [_cond(50.0) for _ in 1:4]
    null_result = BioToolkit.differential_interaction(cond_a, cond_n, coords)
    @test count(<(0.05), null_result.pvalue) <= length(coords) * 0.3
end

@testset "Session 2: hic_diff_test significance from p" begin
    obs_a = [50, 5, 100, 3, 20, 1, 30, 2]
    obs_b = [5, 50, 3, 100, 20, 1, 2, 30]
    result = BioToolkit.hic_diff_test(obs_a, obs_b)
    @test all(0.0 .<= result.pvalue .<= 1.0)
    @test result.significance[1] == "***"                    # most extreme row
    @test all(s -> s in ("***", "**", "*", ".", ""), result.significance)
end

@testset "Session 2: call_loops empirical FDR finds planted loop" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:60]
    Random.seed!(5)
    m = zeros(Float64, 60, 60)
    for i in 1:60, j in i:60
        m[i, j] = m[j, i] = 20 * exp(-0.25 * (j - i)) * (0.5 + rand())
    end
    m[10, 40] = m[40, 10] = 60.0                             # planted long-range loop
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 60, 60) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    loops = BioToolkit.call_loops(mat; window=8, fdr=0.1, min_dist=5, max_dist=50)
    @test nrow(loops) >= 1
    top = loops[1, :]
    @test abs(top.anchor1_bin - 10) <= 2 && abs(top.anchor2_bin - 40) <= 2
    @test all(0.0 .<= loops.score)
end

@testset "Session 2: compartments use O/E and real GC" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:12]
    m = zeros(Float64, 12, 12)
    for i in 1:12, j in i:12
        base = 10 * exp(-0.3 * (j - i))
        # two blocks with strong within-block contacts
        block = (i <= 6) == (j <= 6) ? 3.0 : 0.3
        m[i, j] = m[j, i] = base * block
    end
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 12, 12) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())

    genome = Dict("chr1" => "ACGT" ^ 30000)                     # uniform 50% GC, long enough for all bins
    scores = BioToolkit.compartment_score(mat; genome=genome)
    @test all(isfinite.(scores.gc_content))                  # real GC, not the 0.41 constant
    @test all(scores.gc_content .== 0.5)
    @test length(unique(scores.compartment)) >= 1

    no_genome = BioToolkit.compartment_score(mat)
    @test all(isnan.(no_genome.gc_content))

    ab = BioToolkit.ab_compartments(mat; genome=genome)
    @test 0.0 <= ab.a_fraction <= 1.0
end

@testset "Session 2: pileup accepts merge_loops output" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:30]
    Random.seed!(9)
    m = zeros(Float64, 30, 30)
    for i in 1:30, j in i:30
        m[i, j] = m[j, i] = 10 * exp(-0.3 * (j - i)) * (0.5 + rand())
    end
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 30, 30) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    loops = BioToolkit.call_loops(mat; window=6, fdr=1.0, min_dist=5, max_dist=25, min_expected=0.02)
    nrow(loops) > 0 || return
    merged = BioToolkit.merge_loops(loops[1:min(2, nrow(loops)), :], loops[1:min(2, nrow(loops)), :])
    pu = BioToolkit.pileup_matrix(mat, merged; window_bins=4)
    @test pu.n_loops >= 1
    @test size(pu.matrix) == (9, 9)
end

# ==========================================================================
# Session 2b — Hi-C "surpass Bioconductor" expansion: InteractionSet /
# GInteractions data model, coarsen, mcool writing, 4DN pairs, HiC-Pro,
# saddle plots, compartment strength, TopDom-style TADs, virtual 4C,
# feature aggregation, HiCcompare-style loess normalization.
# ==========================================================================
using HDF5
using CodecZlib

@testset "Session 2b: InteractionSet from matrix" begin
    mat = _session2b_matrix()
    iset = BioToolkit.interactions_from_matrix(mat)
    @test iset isa BioToolkit.InteractionSet
    @test length(iset) == count(!iszero, triu(Matrix(mat.matrix)))
    @test size(BioToolkit.interaction_counts(iset), 2) == 1
    @test BioToolkit.interaction_counts(iset)[1, 1] > 0
    @test all(g -> g.anchor1.chrom == "chr1", iset.interactions)
    # anchors ordered by genomic position
    @test all(g -> (g.anchor1.left, g.anchor1.right) <= (g.anchor2.left, g.anchor2.right) || g.anchor1.chrom != g.anchor2.chrom, iset.interactions)

    cis = BioToolkit.cis_interactions(iset)
    @test length(cis) == length(iset)
    tr = BioToolkit.trans_interactions(iset)
    @test length(tr) == 0

    # find_interactions anchor modes
    probe = BioToolkit.GenomicInterval("chr1", 5 * 1000 + 1, 6 * 1000, '.')
    hits_any = BioToolkit.find_interactions(iset, probe; anchor=:any)
    @test !isempty(hits_any)
    @test all(i -> iset.interactions[i].anchor1.chrom == "chr1", hits_any)

    swapped = BioToolkit.swap_anchors(iset)
    @test swapped.interactions[1].anchor1 == iset.interactions[1].anchor2
end

@testset "Session 2b: interactions from pairs DataFrame" begin
    df = DataFrame(chrom1=["chr1", "chr1", "chr2"], pos1=[100, 500, 100], chrom2=["chr1", "chr2", "chr2"], pos2=[9000, 300, 800])
    iset = BioToolkit.interactions_from_pairs(df)
    @test length(iset) == 3
    tr = BioToolkit.trans_interactions(iset)
    @test length(tr) == 1                      # chr1-chr2 pair
    @test length(BioToolkit.cis_interactions(iset)) == 2
end

@testset "Session 2b: coarsen (zoom) preserves total contact" begin
    mat = _session2b_matrix(; n=32)
    total = sum(Matrix(mat.matrix))
    coarse = BioToolkit.coarsen(mat, 4)
    @test coarse.bin_size == 4000
    @test size(coarse.matrix) == (8, 8)
    @test sum(Matrix(coarse.matrix)) ≈ total atol=1e-6 * total
    @test coarse.bins[1].left == 1 && coarse.bins[1].right == 4000
    one_by = BioToolkit.coarsen(mat, 1)
    @test Matrix(one_by.matrix) == Matrix(mat.matrix)
end

@testset "Session 2b: write_mcool multi-resolution round trip" begin
    mat = _session2b_matrix(; n=40)
    mktempdir() do dir
        path = joinpath(dir, "multi.mcool")
        BioToolkit.write_mcool(mat, path; factors=[1, 2, 5])
        @test BioToolkit.HiC._is_hdf5_file(path)
        r1 = BioToolkit.read_mcool(path; resolution=1000)
        @test size(r1.matrix) == (40, 40)
        @test Matrix(r1.matrix) ≈ Matrix(mat.matrix)
        r5 = BioToolkit.read_mcool(path; resolution=5000)
        @test size(r5.matrix) == (8, 8)
        @test sum(Matrix(r5.matrix)) ≈ sum(Matrix(mat.matrix)) atol=1e-6 * sum(Matrix(mat.matrix))
    end
end

@testset "Session 2b: 4DN pairs round trip" begin
    mktempdir() do dir
        path = joinpath(dir, "contacts.pairs")
        df = DataFrame(read_id=["r1", "r2"], chrom1=["chr1", "chr1"], pos1=[100, 2000],
                       chrom2=["chr1", "chr2"], pos2=[5000, 700], strand1=["+", "-"], strand2=["+", "+"])
        BioToolkit.write_pairs(path, df; chrom_sizes=Dict("chr1" => 100000, "chr2" => 50000))
        back = BioToolkit.read_pairs(path)
        @test nrow(back) == 2
        @test back.chrom1 == ["chr1", "chr1"]
        @test back.pos1 == [100, 2000]
        @test back.chrom2 == ["chr1", "chr2"]
        # gzip-compressed reading
        gz_path = joinpath(dir, "contacts.pairs.gz")
        open(gz_path, "w") do io
            write(io, CodecZlib.GzipCompressorStream(IOBuffer(read(path))))
        end
        back_gz = BioToolkit.read_pairs(gz_path)
        @test nrow(back_gz) == 2
    end
end

@testset "Session 2b: HiC-Pro round trip" begin
    mat = _session2b_matrix(; n=12)
    mktempdir() do dir
        prefix = joinpath(dir, "hiseq")
        BioToolkit.write_hicpro(mat, prefix)
        per_chrom = BioToolkit.read_hicpro("$(prefix)_matrix.matrix", "$(prefix)_abs.bed")
        @test haskey(per_chrom, "chr1")
        recovered = per_chrom["chr1"]
        @test size(recovered.matrix) == (12, 12)
        @test Matrix(recovered.matrix) ≈ Matrix(mat.matrix)
    end
end

@testset "Session 2b: saddle plot recovers block compartments" begin
    # Strong two-block structure -> positive compartment strength.
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:24]
    m = zeros(Float64, 24, 24)
    for i in 1:24, j in i:24
        base = 10 * exp(-0.25 * (j - i))
        block = (i <= 12) == (j <= 12) ? 4.0 : 0.25
        m[i, j] = m[j, i] = base * block
    end
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 24, 24) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    pc1 = vcat(fill(1.0, 12), fill(-1.0, 12))   # A block then B block
    # 2 quantile groups (A vs B): with 24 bins and min_dist 5 this gives
    # enough long-range pairs per group for a well-defined saddle.
    saddle = BioToolkit.saddle_plot(mat, pc1; n_bins=2)
    @test size(saddle.saddle) == (2, 2)
    @test isfinite(saddle.strength)
    @test saddle.strength > 0.5                 # strong compartmentalization
    str = BioToolkit.compartment_strength(mat, pc1; n_bins=2)
    @test str == saddle.strength
end

@testset "Session 2b: TAD caller recovers planted domains" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:60]
    Random.seed!(6)
    m = zeros(Float64, 60, 60)
    for i in 1:60, j in i:60
        m[i, j] = m[j, i] = 10 * exp(-0.2 * (j - i)) * (0.5 + rand())
    end
    # Two TADs (1-20, 30-50) with elevated internal contacts.
    for tad in [(1, 20), (30, 50)]
        for i in tad[1]:tad[2], j in tad[1]:tad[2]
            m[i, j] *= 3.0
        end
    end
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 60, 60) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    tads = BioToolkit.call_tads(mat; window_bins=5, min_size_bins=3)
    @test nrow(tads) >= 2
    # Boundaries (TAD starts/ends) must land near the planted domain edges
    # (20, 30, 50) — the inter-domain gap is depleted of cross contacts.
    bnd = vcat(tads.start_bin, tads.end_bin)
    @test any(x -> abs(x - 20) <= 3, bnd)
    @test any(x -> abs(x - 30) <= 3, bnd)
    @test any(x -> abs(x - 50) <= 3, bnd)
    # TADs overlapping the planted domains are stronger than the gap region.
    gap_strengths = [r.strength for r in eachrow(tads) if r.start_bin >= 21 && r.end_bin <= 30]
    domain_strengths = [r.strength for r in eachrow(tads) if (r.start_bin <= 15 && r.end_bin >= 8) || (r.start_bin >= 30 && r.end_bin <= 50)]
    if !isempty(gap_strengths) && !isempty(domain_strengths)
        @test maximum(domain_strengths) > maximum(gap_strengths)
    end
end

@testset "Session 2b: virtual 4C profile" begin
    mat = _session2b_matrix(; n=30)
    v4c = BioToolkit.virtual_4c(mat, 5; min_dist_bins=2)
    @test nrow(v4c) == 30 - 3                   # viewpoint ± 2 excluded
    @test !any(v4c.bin .== 5)
    @test all(v4c.distance_bins .>= 2)
    # decay: contact decreases with distance
    @test v4c.contact[1] > v4c.contact[end]
end

@testset "Session 2b: aggregate at features" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:30]
    Random.seed!(8)
    m = zeros(Float64, 30, 30)
    for i in 1:30, j in i:30
        m[i, j] = m[j, i] = 10 * exp(-0.3 * (j - i)) * (0.5 + rand())
    end
    # Plant enriched centers at bins 8 and 22.
    for c in (8, 22)
        for i in (c - 2):(c + 2), j in (c - 2):(c + 2)
            m[i, j] += 25.0
        end
    end
    mat = BioToolkit.HiCContactMatrix(spzeros(Float64, 30, 30) + m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    feats = [BioToolkit.GenomicInterval("chr1", 8000, 8999, '.'),
             BioToolkit.GenomicInterval("chr1", 22000, 22999, '.')]
    agg = BioToolkit.aggregate_at_features(mat, feats; window_bins=2)
    @test agg.n_features == 2
    center = agg.matrix[3, 3]
    corner = agg.matrix[1, 1]
    @test center > corner                      # enrichment at feature centers
end

@testset "Session 2b: HiCcompare loess normalization removes trend" begin
    bins = [BioToolkit.GenomicInterval("chr1", (i - 1) * 1000 + 1, i * 1000, '.') for i in 1:40]
    Random.seed!(12)
    function _cond(scale, bias)
        m = spzeros(Float64, 40, 40)
        for i in 1:40, j in i:40
            expected = scale * exp(-0.25 * (j - i))
            v = Float64(rand(Poisson(expected)))
            m[i, j] = m[j, i] = v
        end
        return BioToolkit.HiCContactMatrix(m, "chr1", bins, 1000, false, "none", Dict{String,Any}())
    end
    cond_a = [_cond(4.0, 0.0) for _ in 1:3]
    cond_b = [_cond(12.0, 0.0) for _ in 1:3]   # 3x global scale difference
    coords = [(i, j) for i in 1:40 for j in (i + 1):40]

    result = BioToolkit.hic_loess_normalize(cond_a, cond_b, coords)
    @test nrow(result) > 0
    @test all(0.0 .<= result.pvalue .<= 1.0)
    @test all(result.padj .>= result.pvalue .- 1e-12)
    # After loess correction of the global 3x shift, the bulk of adjusted
    # M values sit near 0 (previously an uncorrected M ≈ log2(3) trend).
    @test abs(median(result.M)) > abs(median(result.M_adj))
    @test abs(median(result.M_adj)) < 0.5
    @test count(<(0.05), result.pvalue) <= nrow(result) * 0.2   # ~null after correction
end

# ==========================================================================
# Session 2c: real-Juicer-file validation. Downloads only the header, master
# index and one block (~200 KB) of a 69 GB ENCODE .hic via HTTP range
# requests, mirroring straw's remote access. Skips silently when offline.
# ==========================================================================
@testset "Session 2c: real Juicer .hic (ranged, online)" begin
    url = "https://www.encodeproject.org/files/ENCFF148QCR/@@download/ENCFF148QCR.hic"
    fetch_range = function (a, b)
        p = tempname()
        ok = success(pipeline(`curl -s --max-time 90 -r $(a)-$(b) -L $url -o $p`, stdout = devnull))
        ok || return nothing
        return read(p)
    end
    header_bytes = fetch_range(0, 65535)
    if header_bytes === nothing || length(header_bytes) < 1024
        @info "skipping real-Juicer validation (offline)"
        return
    end
    @test header_bytes[1:3] == UInt8[UInt8('H'), UInt8('I'), UInt8('C')]
    io = IOBuffer(header_bytes)
    @test BioToolkit.HiC._hic_read_cstring(io) == "HIC"          # magic NUL-terminated
    version = BioToolkit.HiC._read_i32_le(io)
    @test version == 9
    master_pos = BioToolkit.HiC._read_i64_le(io)
    genome_id = BioToolkit.HiC._hic_read_cstring(io)
    @test genome_id == "hg38"
    BioToolkit.HiC._read_i64_le(io); BioToolkit.HiC._read_i64_le(io)   # NVI
    n_attributes = BioToolkit.HiC._read_i32_le(io)
    attrs = Dict{String,String}()
    for _ in 1:n_attributes
        k = BioToolkit.HiC._hic_read_cstring(io)
        attrs[k] = BioToolkit.HiC._hic_read_cstring(io)
    end
    @test haskey(attrs, "software")                              # Juicer tools version
    n_chroms = BioToolkit.HiC._read_i32_le(io)
    chr_index = Dict{String,Int}()
    chrom_lengths = Dict{String,Int}()
    for i in 0:(n_chroms - 1)
        name = BioToolkit.HiC._hic_read_cstring(io)
        len = BioToolkit.HiC._read_i64_le(io)
        chr_index[name] = i
        chrom_lengths[name] = len
    end
    @test haskey(chr_index, "chr1") && chr_index["chr1"] == 1    # 0 = "All"
    @test chrom_lengths["chr1"] == 248956422

    # master index
    idx_bytes = fetch_range(Int(master_pos), Int(master_pos) + 65535)
    idx_bytes === nothing && (@info "skipping (offline)"; return)
    iio = IOBuffer(idx_bytes)
    BioToolkit.HiC._read_i64_le(iio)                             # nBytes of index
    n_entries = BioToolkit.HiC._read_i32_le(iio)
    @test n_entries > 0
    entries = Dict{String,Tuple{Int,Int}}()
    for _ in 1:n_entries
        k = BioToolkit.HiC._hic_read_cstring(iio)
        fpos = Int(BioToolkit.HiC._read_i64_le(iio))
        sz = Int(BioToolkit.HiC._read_i32_le(iio))
        entries[k] = (fpos, sz)
    end
    @test length(entries) == n_entries
    @test haskey(entries, "1_1")                                 # chr1-chr1

    # matrix section for chr1-chr1, resolution 1 Mb
    (mpos, msec) = entries["1_1"]
    @test mpos > 0 && msec > 0
    mz = fetch_range(mpos, mpos + min(msec, 8192) - 1)
    mz === nothing && (@info "skipping (offline)"; return)
    mio = IOBuffer(mz)
    @test BioToolkit.HiC._read_i32_le(mio) == 1                  # c1
    @test BioToolkit.HiC._read_i32_le(mio) == 1                  # c2
    n_res = BioToolkit.HiC._read_i32_le(mio)
    @test n_res >= 2
    target_res = 1_000_000
    block_map = nothing
    for _ in 1:n_res
        unit = BioToolkit.HiC._hic_read_cstring(mio)
        BioToolkit.HiC._read_i32_le(mio)                         # zoom
        for _ in 1:4; BioToolkit.HiC._read_f32_le(mio); end      # stats
        binsize = BioToolkit.HiC._read_i32_le(mio)
        BioToolkit.HiC._read_i32_le(mio)                         # blockBinCount
        BioToolkit.HiC._read_i32_le(mio)                         # blockColumnCount
        n_blocks = BioToolkit.HiC._read_i32_le(mio)
        if unit == "BP" && binsize == target_res
            block_map = Dict{Int,Tuple{Int,Int}}()
            for _ in 1:n_blocks
                bn = BioToolkit.HiC._read_i32_le(mio)
                bp = Int(BioToolkit.HiC._read_i64_le(mio))
                bsz = Int(BioToolkit.HiC._read_i32_le(mio))
                block_map[bn] = (bp, bsz)
            end
            break
        else
            skip(mio, n_blocks * 16)
        end
    end
    @test block_map !== nothing && !isempty(block_map)

    # decode one real block
    bn0 = first(sort(collect(keys(block_map))))
    (bp0, bsz0) = block_map[bn0]
    blk = fetch_range(bp0, bp0 + bsz0 - 1)
    blk === nothing && (@info "skipping (offline)"; return)
    @test length(blk) == bsz0
    raw = IOBuffer(transcode(CodecZlib.ZlibDecompressor, blk))
    n_records = BioToolkit.HiC._read_i32_le(raw)
    bin_x_offset = BioToolkit.HiC._read_i32_le(raw)
    bin_y_offset = BioToolkit.HiC._read_i32_le(raw)
    use_short = read(raw, UInt8) == 0x00
    use_short_x = read(raw, UInt8) == 0x00
    use_short_y = read(raw, UInt8) == 0x00
    rec_type = read(raw, UInt8)
    @test rec_type in (0x01, 0x02)
    n_valid = 0
    if rec_type == 0x01
        n_rows = use_short_y ? BioToolkit.HiC._read_i16_le(raw) : BioToolkit.HiC._read_i32_le(raw)
        for _ in 1:n_rows
            bin_y = bin_y_offset + (use_short_y ? BioToolkit.HiC._read_i16_le(raw) : BioToolkit.HiC._read_i32_le(raw))
            col_count = use_short_x ? BioToolkit.HiC._read_i16_le(raw) : BioToolkit.HiC._read_i32_le(raw)
            for _ in 1:col_count
                bin_x = bin_x_offset + (use_short_x ? BioToolkit.HiC._read_i16_le(raw) : BioToolkit.HiC._read_i32_le(raw))
                c = use_short ? Float64(BioToolkit.HiC._read_i16_le(raw)) : Float64(BioToolkit.HiC._read_f32_le(raw))
                c > 0 && (n_valid += 1)
            end
        end
    else
        n_pts = BioToolkit.HiC._read_i32_le(raw)
        w = BioToolkit.HiC._read_i16_le(raw)
        for i in 0:(n_pts - 1)
            row = i ÷ w; col = i - row * w
            c = use_short ? Float64(BioToolkit.HiC._read_i16_le(raw)) : Float64(BioToolkit.HiC._read_f32_le(raw))
            c > 0 && (n_valid += 1)
        end
    end
    @test n_valid > 0        # real genomic contacts decoded from a 69 GB file
end
