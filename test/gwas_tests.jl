using Test
using SimpleWeightedGraphs
using BioToolkit

@testset "GWAS module" begin
    fixture_prefix = joinpath(@__DIR__, "..", "Examples", "fixtures", "plink_small", "plink_small")
    genotypes = BioToolkit.read_plink(fixture_prefix)
    expected = [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0; 0.0 1.0 1.0]
    @test Matrix(genotypes) == expected
    @test genotypes.bim.snp_id == ["rs1", "rs2", "rs3"]
    @test genotypes.fam.sample_id == ["I1", "I2", "I3", "I4"]

    mktempdir() do dir
        roundtrip_prefix = joinpath(dir, "roundtrip")
        BioToolkit.write_plink(roundtrip_prefix, Matrix(genotypes), genotypes.bim, genotypes.fam)
        roundtrip = BioToolkit.read_plink(roundtrip_prefix)
        @test Matrix(roundtrip) == expected
        @test roundtrip.bim.snp_id == genotypes.bim.snp_id
        @test roundtrip.fam.sample_id == genotypes.fam.sample_id
    end

    phenotype = [0.0, 0.5, 1.0, 1.5]

    linear = BioToolkit.gwas_linear_scan(genotypes, phenotype)
    @test linear.method == "linear_scan"
    @test length(linear.snp_ids) == 3
    @test all(isfinite, linear.pvalue)

    mixed = BioToolkit.gwas_lmm_scan(genotypes, phenotype)
    @test mixed.method == "linear_scan"
    @test length(mixed.snp_ids) == 3

    clumped = BioToolkit.ld_clumping(linear, genotypes; p_threshold=1.0, r2_threshold=0.1, window=500)
    @test length(clumped.snp_ids) <= length(linear.snp_ids)

    prs = BioToolkit.calculate_prs(genotypes, linear)
    @test length(prs) == size(genotypes, 1)

    prs_cv = BioToolkit.prs_cross_validation(genotypes, phenotype, linear)
    @test haskey(prs_cv.scores, (1e-5, 50_000))

    meta = BioToolkit.meta_analyze([linear, linear])
    @test meta.study_count == 2
    @test length(meta.snp_ids) == length(linear.snp_ids)

    peakset = BioToolkit.PeakSet([
        BioToolkit.Peak("peak1", "chr1", 50, 150, 100, 2.0, 1e-6, 1e-4),
    ])
    overlaps = BioToolkit.overlap_gwas_peaks(linear, peakset; flank=75, pvalue_threshold=1.0)
    @test "snp_id" in names(overlaps)

    network = BioToolkit.GeneNetwork(
        SimpleWeightedGraph(3),
        Dict("GENE1" => 1, "GENE2" => 2, "GENE3" => 3),
        ["GENE1", "GENE2", "GENE3"],
        [0.0, 0.0, 0.0],
        [1, 1, 2],
    )
    enrichment = BioToolkit.gene_based_test(linear, network; pvalue_threshold=1.0)
    @test enrichment isa Vector

    manhattan_points, _ = BioToolkit.BioPlotting.manhattan_data(linear)
    qq_points, _ = BioToolkit.BioPlotting.qq_data(linear)
    forest = BioToolkit.BioPlotting.gwas_forest_plot(meta)
    manhattan = BioToolkit.BioPlotting.manhattan_plot(linear)
    qq = BioToolkit.BioPlotting.qq_plot(linear)
    @test !isempty(manhattan_points)
    @test !isempty(qq_points)
    @test forest.figure !== nothing
    @test manhattan.figure.subplots[1].attr[:title] == "Manhattan plot"
    @test manhattan.figure.subplots[1].attr[:annotations] |> length == 3
    @test length(manhattan.figure.series_list) == 2
    @test qq.figure.subplots[1].attr[:title] == "QQ plot"
    @test qq.figure.subplots[1].attr[:aspect_ratio] == :equal
    @test length(qq.figure.series_list) == 3
end