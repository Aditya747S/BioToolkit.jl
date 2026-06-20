using Test
using Random
using DataFrames
using Statistics
using BioToolkit

@testset "Bioconductor compatibility wrappers" begin
    @testset "Single-cell wrappers" begin
        counts = [
            12 10 2 1 0 0;
            8 9 1 1 0 0;
            1 0 8 9 2 1;
            0 0 2 1 9 10
        ]
        gene_ids = ["MT-CO1", "MT-ND1", "GeneA", "GeneB"]
        cell_ids = ["c1", "c2", "c3", "c4", "c5", "c6"]
        groups = ["A", "A", "B", "B", "C", "C"]

        sce = BioToolkit.SingleCellExperiment(counts, gene_ids, cell_ids)

        qc = BioToolkit.scater_qc_metrics(sce)
        @test nrow(qc) == 6
        @test all(name in names(qc) for name in ["cell_id", "total_counts", "detected_features", "mito_counts", "pct_mito"])
        @test all(isfinite, qc.pct_mito)

        sf = BioToolkit.scran_pool_size_factors(sce; pool_size=3, overlap=2)
        @test length(sf) == 6
        @test all(isfinite, sf)
        @test all(>(0.0), sf)
        @test isapprox(median(sf), 1.0; atol=1e-6)

        pb = BioToolkit.scuttle_aggregate_counts(sce, groups)
        @test pb.sample_ids == ["A", "B", "C"]
        @test size(pb.counts) == (4, 3)

        annotation = BioToolkit.singler_annotate(sce; labels=groups)
        @test length(annotation.cell_labels) == length(cell_ids)

        projection = BioToolkit.scmap_project(sce, sce; labels=groups, n_components=2)
        @test size(projection.embedding, 1) == length(cell_ids)
        @test size(projection.embedding, 2) == 2
        @test length(projection.predicted_label) == length(cell_ids)

        sim = BioToolkit.splatter_simulate_counts(n_genes=80, n_cells=40, n_groups=2, seed=11)
        @test size(sim.experiment.counts) == (80, 40)
        @test nrow(sim.truth) == 80
    end

    @testset "Epigenetics wrappers" begin
        fragments = Dict(
            "s1" => [(chrom="chr1", left=100, right=180), (chrom="chr1", left=130, right=210), (chrom="chr1", left=150, right=220)],
            "s2" => [(chrom="chr1", left=110, right=190), (chrom="chr1", left=145, right=230)],
            "s3" => [(chrom="chr1", left=500, right=580), (chrom="chr1", left=520, right=600), (chrom="chr1", left=540, right=620)],
            "s4" => [(chrom="chr1", left=510, right=590), (chrom="chr1", left=550, right=630)]
        )
        design = [:control, :control, :treated, :treated]
        csaw = BioToolkit.csaw_window_differential_binding(fragments, design; window_size=120, step_size=60, min_count=1)
        @test nrow(csaw) > 0
        @test all(name in names(csaw) for name in ["gene_id", "chromosome", "start", "stop", "pvalue", "padj"])

        beta = [
            0.85 0.88 0.41 0.43;
            0.20 0.25 0.65 0.70;
            0.50 0.55 0.49 0.46
        ]
        dmp = BioToolkit.minfi_dmp(beta, [:A, :A, :B, :B]; feature_ids=["cg1", "cg2", "cg3"])
        @test nrow(dmp) == 3
        @test all(isfinite, dmp.delta_beta)
        @test all(x -> ismissing(x) || (0.0 <= x <= 1.0), dmp.padj)

        calls = BioToolkit.MethylationCall[
            BioToolkit.MethylationCall("chr1", 100, "A1", true),
            BioToolkit.MethylationCall("chr1", 100, "A2", true),
            BioToolkit.MethylationCall("chr1", 100, "B1", false),
            BioToolkit.MethylationCall("chr1", 100, "B2", false),
            BioToolkit.MethylationCall("chr1", 150, "A1", true),
            BioToolkit.MethylationCall("chr1", 150, "A2", false),
            BioToolkit.MethylationCall("chr1", 150, "B1", false),
            BioToolkit.MethylationCall("chr1", 150, "B2", false)
        ]
        dmr = BioToolkit.bsseq_dmr(calls, [:A, :A, :B, :B]; bin_width=100, min_total=1)
        @test isa(dmr, DataFrame)
        @test all(name in names(dmr) for name in ["region_id", "chrom", "start", "stop", "pvalue", "padj"])
    end

    @testset "Proteomics, metabolomics, multi-omics" begin
        rt = collect(1.0:1.0:120.0)
        mz1 = fill(500.0, length(rt))
        mz2 = fill(500.005, length(rt))
        intensity1 = [exp(-((t - 35.0)^2) / 50.0) + 0.8exp(-((t - 85.0)^2) / 70.0) for t in rt]
        intensity2 = [1.1exp(-((t - 37.0)^2) / 55.0) + 0.7exp(-((t - 83.0)^2) / 65.0) for t in rt]

        spectra = [
            BioToolkit.Spectrum(1.0, mz1, intensity1),
            BioToolkit.Spectrum(2.0, mz2, intensity2)
        ]
        ms = BioToolkit.MassSpecExperiment(spectra)
        xcms = BioToolkit.xcms_peak_workflow(ms; scales=[1.0, 2.0], threshold=0.4)
        @test nrow(xcms.peak_table) >= 1
        @test size(xcms.alignment_cost, 1) == 2

        matrix = [
            10.0 11.0 14.0 15.0;
            8.0 8.5 12.0 11.5;
            5.0 4.5 5.2 5.0
        ]
        da = BioToolkit.lipidr_differential_abundance(matrix, [:ctrl, :ctrl, :case, :case])
        @test size(da.coefficients, 1) == size(matrix, 1)

        Random.seed!(5)
        assay1 = randn(12, 25)
        assay2 = randn(12, 15)
        mix = BioToolkit.mixomics_factor_analysis([assay1, assay2]; n_components=3)
        mofa = BioToolkit.mofa2_factor_analysis([assay1, assay2]; n_factors=3)
        @test size(mix.factors, 2) == 3
        @test size(mofa.factors, 2) == 3
    end

    @testset "Variant and visualization wrappers" begin
        variants = [
            (chrom="chr1", pos=100, ref="A", alt="G", qual=40.0, filter="PASS"),
            (chrom="chr1", pos=120, ref="AT", alt="A", qual=60.0, filter="PASS"),
            (chrom="chr1", pos=130, ref="C", alt="T", qual=10.0, filter="q10")
        ]

        filtered = BioToolkit.varianttools_filter_variants(variants; min_qual=20.0, pass_only=true, snv_only=true)
        @test length(filtered) == 1
        @test filtered[1].pos == 100

        heat = BioToolkit.complexheatmap_payload(randn(6, 4); row_labels=["g$(i)" for i in 1:6], column_labels=["s$(i)" for i in 1:4])
        @test size(heat.matrix) == (6, 4)
        @test length(heat.row_order) == 6

        track = BioToolkit.gviz_track_table(["chr1", "chr2"], [100, 200], [150, 280]; scores=[3.2, 1.5], labels=["peak1", "peak2"], track="coverage")
        @test nrow(track) == 2
        @test all(name in names(track) for name in ["track", "chromosome", "start", "stop", "score", "label"])
    end

    @testset "Parallel wrapper" begin
        values = BioToolkit.biocparallel_map(x -> x^2, 1:12; nworkers=3)
        @test values == [i^2 for i in 1:12]
    end
end
