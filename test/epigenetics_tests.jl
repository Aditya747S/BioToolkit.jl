using SparseArrays
using DataFrames

@testset "Epigenetics" begin
    fragments = [
        BioToolkit.GenomicInterval("chr1", 1, 5),
        BioToolkit.GenomicInterval("chr1", 4, 8),
        BioToolkit.GenomicInterval("chr2", 2, 4),
    ]

    coverage = BioToolkit.calculate_coverage(fragments; chrom_lengths=Dict("chr1" => 10, "chr2" => 6))
    @test BioToolkit.coverage_depth(coverage["chr1"], 1) == 1
    @test BioToolkit.coverage_depth(coverage["chr1"], 4) == 2
    @test BioToolkit.coverage_depth(coverage["chr1"], 9) == 0
    @test !isempty(BioToolkit.coverage_segments(coverage["chr1"]; chrom="chr1"))

    peaks = BioToolkit.call_peaks(coverage; pvalue_threshold=1.0, min_depth=1)
    @test !isempty(peaks)
    @test all(peak -> peak.chrom in ("chr1", "chr2"), peaks.peaks)

    merged_coverage = BioToolkit.calculate_coverage([
        BioToolkit.GenomicInterval("chr1", 1, 3),
        BioToolkit.GenomicInterval("chr1", 5, 7),
    ]; chrom_lengths=Dict("chr1" => 10))
    merged_peaks = BioToolkit.call_peaks(merged_coverage; pvalue_threshold=1.0, min_depth=1, merge_gap=1)
    @test length(merged_peaks.peaks) == 1
    @test merged_peaks.peaks[1].left == 1
    @test merged_peaks.peaks[1].right == 7

    fragments_by_sample = Dict(
        "s1" => fragments[1:2],
        "s2" => fragments[3:3],
    )
    count_matrix = BioToolkit.count_overlaps(fragments_by_sample, peaks)
    @test size(count_matrix.counts, 2) == 2
    @test !isempty(count_matrix.gene_ids)

    binding = BioToolkit.differential_binding(fragments_by_sample, peaks, [:control, :treated]; min_total=0, shrink=false)
    @test length(binding) == size(count_matrix.counts, 1)

    support = BioToolkit.summarize_peak_support(fragments_by_sample, peaks)
    @test !isempty(support)
    @test any(item -> item.fragment_count > 0, support)
    @test all(item -> isfinite(item.enrichment), support)

    epigenome = BioToolkit.Epigenome(
        fragments,
        coverage,
        count_matrix.counts,
        DataFrame(sample_id=["s1", "s2"], condition=["control", "treated"]),
    )
    @test length(epigenome.intervals) == 3
    @test haskey(epigenome.coverage, "chr1")

    gc_coverage = BioToolkit.SparseCoverageVector([1, 3, 5, 7, 9], [10, 8, 2, 1, 0], 8)
    gc_result = BioToolkit.normalize_gc_bias("GGGGAAAA" , gc_coverage; window_size=2)
    @test length(gc_result.gc) == 4
    @test gc_result.corrected_coverage isa BioToolkit.SparseCoverageVector
    @test gc_result.corrected_coverage.depths[1] < gc_result.observed[1]
    @test gc_result.corrected_coverage.depths[3] > round(Int, gc_result.observed[3])

    calls = [
        BioToolkit.MethylationCall("chr1", 1, "s1", true),
        BioToolkit.MethylationCall("chr1", 2, "s1", false),
        BioToolkit.MethylationCall("chr1", 1, "s2", true),
        BioToolkit.MethylationCall("chr1", 2, "s2", true),
        BioToolkit.MethylationCall("chr1", 1, "s3", false),
        BioToolkit.MethylationCall("chr1", 2, "s3", false),
        BioToolkit.MethylationCall("chr1", 1, "s4", true),
        BioToolkit.MethylationCall("chr1", 2, "s4", false),
    ]
    methylation = BioToolkit.bin_methylation(calls; bin_width=10, sample_metadata=DataFrame(sample_id=["s1", "s2", "s3", "s4"]))
    @test size(methylation.methylated, 2) == 4
    @test length(methylation.sample_ids) == 4
    methylation_results = BioToolkit.differential_methylation(methylation, [:control, :control, :treated, :treated]; min_total=0)
    @test length(methylation_results) == size(methylation.total, 1)
    @test all(result -> 0.0 <= result.group1_mean <= 1.0 && 0.0 <= result.group2_mean <= 1.0, methylation_results)
    @test all(result -> isfinite(result.stat) && isfinite(result.pvalue) && isfinite(result.padj), methylation_results)

    base = BioToolkit.SingleCellExperiment([1 0 1 0; 0 1 0 1; 1 1 0 0], ["g1", "g2", "g3"], ["c1", "c2", "c3", "c4"])
    chromatin = BioToolkit.SingleCellChromatinExperiment(base, sparse([1 0 1 0; 0 1 1 0; 1 1 0 1]), ["p1", "p2", "p3"], Dict{String,Matrix{Float64}}(), Dict{String,Any}())
    tfidf = BioToolkit.tfidf(chromatin)
    @test size(tfidf) == (3, 4)
    lsi = BioToolkit.run_lsi(chromatin; n_components=2)
    @test size(lsi, 2) == 2
    @test size(BioToolkit.rsvd(tfidf; rank=2)[1], 2) == 2

    peak_intervals = [
        BioToolkit.GenomicInterval("chr1", 1, 5),
        BioToolkit.GenomicInterval("chr1", 4, 8),
        BioToolkit.GenomicInterval("chr1", 7, 9),
    ]
    gene_intervals = [BioToolkit.GenomicInterval("chr1", 1, 9, '+')]
    activity = BioToolkit.gene_activity_score(chromatin, gene_intervals, peak_intervals)
    @test size(activity) == (1, 4)

    coaccessibility = BioToolkit.calculate_coaccessibility(chromatin; min_correlation=-1.0)
    @test all(edge -> edge isa BioToolkit.CoaccessibilityEdge, coaccessibility)

    motif_deviation = BioToolkit.compute_motif_deviations(chromatin, Dict("motif1" => [1, 2]))
    @test haskey(motif_deviation, "motif1")

    footprint = BioToolkit.detect_footprints([0.2, 0.4, 0.1, 0.05, 0.3, 0.45, 0.2], [4]; flank=2, center=1)
    @test length(footprint) == 1

    contact = BioToolkit.ContactMatrix(sparse([0 4 1 0; 4 0 2 1; 1 2 0 3; 0 1 3 0]), peak_intervals, 1000, "chr1")
    di = BioToolkit.directionality_index(contact; window=1)
    @test length(di) == 4
    tads = BioToolkit.detect_tads(contact; window=1, threshold=0.1)
    @test all(tad -> tad isa BioToolkit.TadResult, tads)
    @test !isempty(tads)
end