using Test
using DataFrames

function _make_browser_bam(dir)
    bam_path = joinpath(dir, "browser.bam")
    header = BioToolkit.BamHeader([BioToolkit.BamReference("chr1", 1_000)])
    records = [
        BioToolkit.BamRecord("read1", "chr1", 99, [BioToolkit.BamCigarOp(10, 'M')], "ACGTACGTAA"; tags=Dict("MD" => "5A4")),
        BioToolkit.BamRecord("read2", "chr1", 119, [BioToolkit.BamCigarOp(8, 'M')], "TTTTGGGG"; mate_refname="chr1", mate_pos=159, template_length=40),
        BioToolkit.BamRecord("read3", "chr1", 300, [BioToolkit.BamCigarOp(6, 'M')], "CCCCCC"),
    ]
    BioToolkit.write_bam(bam_path, BioToolkit.BamFile(records; header=header))
    return bam_path
end

@testset "Genome browser core" begin
    viewport = BioToolkit.GenomeViewport("chr1", 100, 499; pixel_width=12)
    @test viewport.chrom == "chr1"
    @test first(viewport.range) == 100
    @test last(viewport.range) == 499
    @test viewport.pixel_width == 12
    @test viewport.resolution ≈ (499 - 100) / 12
    @test BioToolkit.genome_lod(viewport) == :gene

    chromosome_view = BioToolkit.GenomeViewport("chr1", 1, 100_000; pixel_width=5)
    @test BioToolkit.genome_lod(chromosome_view) == :chromosome

    basepair_view = BioToolkit.GenomeViewport("chr1", 1, 180; pixel_width=20)
    @test BioToolkit.genome_lod(basepair_view) == :basepair
end

@testset "Gene track packing" begin
    genes = [
        BioToolkit.GenomicInterval("chr1", 100, 160, '+', Dict("name" => "geneA", "exons" => [(100, 120), (140, 160)])),
        BioToolkit.GenomicInterval("chr1", 120, 190, '-', Dict("name" => "geneB")),
        BioToolkit.GenomicInterval("chr1", 220, 250, '+', Dict("name" => "geneC")),
    ]

    track = BioToolkit.GeneTrack(genes; label_field=:name, max_labels=2)
    chromosome_view = BioToolkit.GenomeViewport("chr1", 1, 100_000; pixel_width=5)
    chromosome_plan = BioToolkit.render_track(track, chromosome_view)
    @test chromosome_plan.lod == :chromosome
    @test !isempty(chromosome_plan.density)

    gene_view = BioToolkit.GenomeViewport("chr1", 90, 260; pixel_width=5)
    gene_plan = BioToolkit.render_track(track, gene_view)
    @test gene_plan.lod == :gene
    @test gene_plan.row_count == 2
    @test [placement.row for placement in gene_plan.placements] == [1, 2, 1]
    @test gene_plan.labels == ["geneA", "geneB"]

    basepair_view = BioToolkit.GenomeViewport("chr1", 90, 160; pixel_width=20)
    basepair_plan = BioToolkit.render_track(track, basepair_view)
    @test basepair_plan.lod == :basepair
    @test any(segment -> segment.kind == :exon, basepair_plan.placements[1].segments)
end

@testset "Coverage track aggregation" begin
    sparse = BioToolkit.SparseCoverageVector([1, 3, 5], [2, 4, 6], 8)
    track = BioToolkit.CoverageTrack(sparse; chrom="chr1", start=100, window_size=2)
    viewport = BioToolkit.GenomeViewport("chr1", 100, 500; pixel_width=10)
    plan = BioToolkit.render_track(track, viewport)
    @test plan.lod == :gene
    @test length(plan.bins) == 4
    @test plan.bins[1].mean_value == 1.0
    @test plan.bins[2].max_value == 4.0

    table = DataFrame(chrom=["chr1", "chr1", "chr1"], start=[100, 110, 120], stop=[105, 115, 125])
    table_track = BioToolkit.CoverageTrack(table; sorted=true, window_size=10)
    table_plan = BioToolkit.render_track(table_track, BioToolkit.GenomeViewport("chr1", 100, 129; pixel_width=3))
    @test !isempty(table_plan.bins)
    @test all(bin -> bin.min_value == bin.max_value == bin.mean_value, table_plan.bins)
end

@testset "Alignment track rendering" begin
    mktempdir() do dir
        bam_path = _make_browser_bam(dir)
        track = BioToolkit.AlignmentTrack(bam_path; show_mismatches=true, max_reads=10, rng_seed=7)
        viewport = BioToolkit.GenomeViewport("chr1", 100, 170; pixel_width=20)
        plan = BioToolkit.render_track(track, viewport)
        @test plan.lod == :basepair
        @test length(plan.reads) == 2
        @test any(read -> !isempty(read.mismatches), plan.reads)
        @test length(plan.connectors) == 1
        @test plan.connectors[1].qname == "read2"

        sampled_track = BioToolkit.AlignmentTrack(bam_path; show_mismatches=true, max_reads=1, rng_seed=1)
        sampled_plan = BioToolkit.render_track(sampled_track, viewport)
        @test sampled_plan.downsampled == true
        @test length(sampled_plan.reads) == 1
    end
end

@testset "Browser composition" begin
    genes = [BioToolkit.GenomicInterval("chr1", 100, 160, '+', Dict("name" => "geneA"))]
    track1 = BioToolkit.GeneTrack(genes)
    track2 = BioToolkit.CoverageTrack([0.0, 1.0, 2.0, 3.0]; chrom="chr1", start=100)
    browser = BioToolkit.GenomeBrowser([track1, track2], BioToolkit.GenomeViewport("chr1", 100, 160; pixel_width=40); title="Test Browser")
    plans = BioToolkit.render_browser(browser)
    @test length(plans) == 2
    @test plans[1].plan isa BioToolkit.GeneRenderPlan
    @test plans[2].plan isa BioToolkit.CoverageRenderPlan
    @test occursin("GenomeBrowser", sprint(show, browser))
end

@testset "Makie smoke test" begin
    if Base.find_package("Makie") !== nothing && Base.find_package("CairoMakie") !== nothing
        @eval using Makie
        @eval using CairoMakie

        figure = Makie.Figure()
        axis = Makie.Axis(figure[1, 1])
        Makie.lines!(axis, [1, 2, 3], [1, 4, 9])
        @test figure isa Makie.Figure

        mktempdir() do dir
            output_path = joinpath(dir, "browser_smoke.png")
            CairoMakie.save(output_path, figure)
            @test isfile(output_path)
            @test filesize(output_path) > 0
        end
    else
        @test true
    end
end