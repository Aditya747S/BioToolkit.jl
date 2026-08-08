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

    @testset "ATAC-seq" begin
        atac_fragments = [
            BioToolkit.GenomicInterval("chr1", 1, 100),
            BioToolkit.GenomicInterval("chr1", 150, 250),
            BioToolkit.GenomicInterval("chr1", 300, 400),
            BioToolkit.GenomicInterval("chr1", 450, 550),
            BioToolkit.GenomicInterval("chr1", 10000, 10100),
            BioToolkit.GenomicInterval("chr1", 10050, 10150),
            BioToolkit.GenomicInterval("chr1", 200, 300),
            BioToolkit.GenomicInterval("chr1", 350, 450),
        ]
        barcodes = ["cell1", "cell2", "cell1", "cell3", "cell4", "cell4", "cell1", "cell2"]
        sample_ids = ["sample1", "sample1", "sample1", "sample1", "sample1", "sample1", "sample1", "sample1"]

        exp = BioToolkit.atac_experiment(atac_fragments; barcodes=barcodes, sample_ids=sample_ids)
        @test length(exp) == 8
        @test !isempty(exp)
        @test exp.sample_ids == ["sample1"]

        frags = BioToolkit.ATACFragment[]
        for i in 1:20
            push!(frags, BioToolkit.ATACFragment("chr1", i * 100, i * 100 + 100, 100, '+', "cell_$i", "sample1"))
        end
        for i in 1:10
            push!(frags, BioToolkit.ATACFragment("chr1", 147 + i * 200, 147 + i * 200 + 147, 147, '+', "cell_nuc_$i", "sample1"))
        end
        atac_exp = BioToolkit.ATACExperiment(frags, Dict("chr1" => collect(1:length(frags))), ["sample1"], Dict{String,Any}())

        frag_dist = BioToolkit.fragment_size_distribution(atac_exp)
        @test frag_dist isa BioToolkit.FragmentSizeDistribution
        @test frag_dist.nucleosome_free > 0.0
        @test frag_dist.mono_nucleosome > 0.0
        @test frag_dist.mean_size > 0.0
        @test frag_dist.median_size > 0.0
        @test sum(frag_dist.counts) == length(frags)

        nuc_metrics = BioToolkit.nucleosome_metrics(atac_exp)
        @test nuc_metrics isa BioToolkit.NucleosomeMetricsResult
        @test nuc_metrics.nucleosome_signal >= 0.0
        @test nuc_metrics.mononucleosome_fraction >= 0.0
        @test nuc_metrics.fraction_small >= 0.0

        tss_sites = [
            BioToolkit.GenomicInterval("chr1", 1000, 1001, '+'),
            BioToolkit.GenomicInterval("chr1", 5000, 5001, '+'),
            BioToolkit.GenomicInterval("chr1", 10000, 10001, '+'),
        ]

        tss_result = BioToolkit.tss_enrichment(atac_exp, tss_sites; flank=500, bin_width=50)
        @test haskey(tss_result, :positions)
        @test haskey(tss_result, :profile)
        @test haskey(tss_result, :enrichment)
        @test length(tss_result.positions) == length(tss_result.profile)
        @test tss_result.enrichment >= 0.0

        peaks_for_frip = BioToolkit.PeakSet([
            BioToolkit.Peak("peak1", "chr1", 100, 200, 150, 10.0, 0.01, 0.05),
            BioToolkit.Peak("peak2", "chr1", 10000, 10100, 10050, 8.0, 0.02, 0.08),
        ])

        frip = BioToolkit.frip_score(atac_exp, peaks_for_frip)
        @test frip >= 0.0
        @test frip <= 1.0

        bias_result = BioToolkit.insertion_bias_correction(atac_exp)
        @test haskey(bias_result, :correction_factors)
        @test bias_result.total_insertions > 0

        atac_peaks = BioToolkit.atac_peak_calling(atac_exp; pvalue_threshold=1.0, min_depth=1)
        @test atac_peaks isa BioToolkit.PeakSet

        qc_result = BioToolkit.atac_qc_report(atac_exp, peaks_for_frip, tss_sites)
        @test qc_result isa BioToolkit.ATACQCResult
        @test qc_result.total_fragments > 0
        @test qc_result.frip >= 0.0
        @test qc_result.tss_enrichment >= 0.0
        @test qc_result.nucleosome_signal >= 0.0
        @test 0.0 <= qc_result.fraction_nucleosome_free <= 1.0
        @test 0.0 <= qc_result.fraction_mono_nucleosome <= 1.0
    end
end