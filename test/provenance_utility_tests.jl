using Test
using DataFrames
using Tables
using Graphs

@testset "Provenance utilities" begin
    @testset "Protein" begin
        ctx = BioToolkit.ProvenanceContext()
        seq = BioToolkit.AASeq("ACDEFGHIK")

        mass = BioToolkit.protein_mass(seq; prov_ctx=ctx)
        ec = BioToolkit.extinction_coefficient(seq; prov_ctx=ctx)
        pi = BioToolkit.isoelectric_point(seq; prov_ctx=ctx)
        gravy = BioToolkit.gravy(seq; prov_ctx=ctx)
        ai = BioToolkit.aliphatic_index(seq; prov_ctx=ctx)
        summary = BioToolkit.protparam(seq; prov_ctx=ctx)

        @test mass > 0
        @test ec >= 0
        @test 0.0 <= pi <= 14.0
        @test isfinite(gravy)
        @test ai >= 0
        @test summary.length == length(seq)

        @test _node_by_operation(ctx, "protein_mass").parameters["type"] == "monoisotopic"
        @test haskey(_node_by_operation(ctx, "protparam").parameters, "molecular_weight_mono")
    end

    @testset "Quality" begin
        ctx = BioToolkit.ProvenanceContext()
        fastq = BioToolkit.FastqRecord(BioToolkit.DNASeq("ACGTACGT"), "IIIIIIII"; identifier="r1")
        lite = BioToolkit.SeqRecordLite(BioToolkit.DNASeq("ACGTACGT"); identifier="r2", letter_annotations=Dict(:quality => "IIIIIIII"))

        scores = BioToolkit.phred_scores("IIII"; prov_ctx=ctx)
        roundtrip = BioToolkit.phred_string(scores; prov_ctx=ctx)
        mean_q = BioToolkit.mean_quality("IIII"; prov_ctx=ctx)
        ok = BioToolkit.quality_filter(fastq; prov_ctx=ctx)
        trimmed = BioToolkit.trim_low_quality(fastq; prov_ctx=ctx)
        processed = BioToolkit.process_sequencing_record(lite; prov_ctx=ctx)
        batch = BioToolkit.sequencing_pipeline([fastq, fastq]; prov_ctx=ctx)

        @test scores == [40, 40, 40, 40]
        @test roundtrip == "IIII"
        @test mean_q == 40.0
        @test ok == true
        @test trimmed !== nothing
        @test processed !== nothing
        @test length(batch) == 2

        @test _node_by_operation(ctx, "phred_scores").parameters["score_count"] == 4
        @test _node_by_operation(ctx, "quality_filter").parameters["result"] == true
        @test _node_by_operation(ctx, "sequencing_pipeline").parameters["output_count"] == 2
    end

    @testset "Query" begin
        ctx = BioToolkit.ProvenanceContext()
        table = (
            chrom=["chr1", "chr1", "chr2"],
            pos=[10, 120, 30],
            id=["a", "b", "c"],
            ref=["A", "C", "G"],
            alt=["T", "G", "A"],
            qual=Union{Missing,Float32}[missing, Float32(10), Float32(20)],
        )

        subset = BioToolkit.filter_region(table, "chr1", 0, 100; prov_ctx=ctx)
        bins = BioToolkit.bin_positions(Int32[1, 2, 11, 12], 10; prov_ctx=ctx)
        hist = BioToolkit.coverage_histogram(table, "chr1", 10; prov_ctx=ctx)
        cov = BioToolkit.window_coverage(Int32[1, 5, 9], Int32[3, 8, 12], 4; prov_ctx=ctx)

        tmp = mktempdir()
        out1 = joinpath(tmp, "cov.bw")
        out2 = joinpath(tmp, "cov_map.bw")
        BioToolkit.write_bigwig(Float64[0, 1, 1, 2], out1; prov_ctx=ctx)
        BioToolkit.write_bigwig(Dict("chr1" => Float64[0, 1, 2]), out2; prov_ctx=ctx)

        @test subset.pos == [10]
        @test bins[0] == 2
        @test hist[1] == 1
        @test !isempty(cov)
        @test isfile(out1)
        @test isfile(out2)

        @test _node_by_operation(ctx, "filter_region").parameters["chrom"] == "chr1"
        @test _node_by_operation(ctx, "bin_positions").parameters["bin_count"] >= 1
        write_nodes = [node for node in values(ctx.nodes) if node.operation == "write_bigwig"]
        @test any(node.parameters["output"] == out1 for node in write_nodes)
        @test any(node.parameters["output"] == out2 for node in write_nodes)
    end

    @testset "Browser" begin
        ctx = BioToolkit.ProvenanceContext()
        viewport = BioToolkit.GenomeViewport("chr1", 100, 499; pixel_width=12, prov_ctx=ctx)
        genes = [BioToolkit.GenomicInterval("chr1", 120, 180, '+', Dict("name" => "geneA"))]
        gene_track = BioToolkit.GeneTrack(genes; prov_ctx=ctx)
        cov_track = BioToolkit.CoverageTrack([0.0, 1.0, 2.0, 3.0]; chrom="chr1", start=100, prov_ctx=ctx)
        browser = BioToolkit.GenomeBrowser([gene_track, cov_track], viewport; title="Browser", gap=8)
        plans = BioToolkit.render_browser(browser; prov_ctx=ctx)

        @test viewport.pixel_width == 12
        @test length(plans) == 2
        @test plans[1].plan isa BioToolkit.GeneRenderPlan
        @test _node_by_operation(ctx, "GenomeViewport").parameters["pixel_width"] == 12
        @test _node_by_operation(ctx, "render_browser").parameters["track_count"] == 2
    end

    @testset "Record and schema" begin
        ctx = BioToolkit.ProvenanceContext()
        dna = BioToolkit.DNASeq("ACGTACGT")

        seq_record = BioToolkit.SeqRecord(dna; identifier="seq1", prov_ctx=ctx)
        fastq_record = BioToolkit.FastqRecord(dna, "IIIIIIII"; identifier="fq1", prov_ctx=ctx)
        lite_record = BioToolkit.SeqRecordLite(dna; identifier="lite1", letter_annotations=Dict(:quality => "IIIIIIII"), prov_ctx=ctx)
        variant = BioToolkit.VariantTextRecord("chr1", 10, "var1", "A", "T", 12.0; prov_ctx=ctx)
        document = BioToolkit.VcfDocument([variant]; prov_ctx=ctx)
        schema = BioToolkit.arrow_schema(BioToolkit.VariantEvent; prov_ctx=ctx)
        compact = BioToolkit.compact_variant_event(variant; prov_ctx=ctx)

        @test length(seq_record) == 8
        @test length(fastq_record) == 8
        @test length(lite_record) == 8
        @test length(document) == 1
        @test schema isa Tables.Schema
        @test compact.hasqual == true

        @test _node_by_operation(ctx, "SeqRecord").parameters["identifier"] == "seq1"
        @test _node_by_operation(ctx, "FastqRecord").parameters["length"] == 8
        @test _node_by_operation(ctx, "VcfDocument").parameters["record_count"] == 1
        @test _node_by_operation(ctx, "arrow_schema").parameters["field_count"] == length(fieldnames(BioToolkit.VariantEvent))
    end

    @testset "Proteomics" begin
        ctx = BioToolkit.ProvenanceContext()
        matrix = [1.0 2.0; 3.0 4.0]
        groups = ["A", "B"]

        imputed = BioToolkit.qrilc_impute(matrix; prov_ctx=ctx)
        da = BioToolkit.differential_abundance(matrix, groups; prov_ctx=ctx)
        dia = BioToolkit.dia_like_quantification(matrix, groups; prov_ctx=ctx)
        motifs = BioToolkit.glycoproteomics_motif_table(["NVT", "AAA"]; prov_ctx=ctx)
        peptides = BioToolkit.protein_inference_top3(DataFrame(protein=["P1", "P1", "P2"], intensity=[1.0, 2.0, 3.0]); prov_ctx=ctx)
        enrichment = BioToolkit.ptm_site_enrichment(["S1", "S2"], ["S1", "S1", "S2"]; prov_ctx=ctx)
        graph, labels, peaks = BioToolkit.build_spectrum_graph([100.0, 171.0, 228.0]; prov_ctx=ctx)

        @test size(imputed) == size(matrix)
        @test da.coefficients isa Matrix{Float64}
        @test nrow(dia) == 2
        @test nrow(motifs) == 2
        @test nrow(peptides) == 2
        @test nrow(enrichment) >= 1
        @test graph isa SimpleDiGraph
        @test length(labels) >= 1
        @test length(peaks) == 3

        @test _node_by_operation(ctx, "qrilc_impute").parameters["row_count"] == 2
        @test _node_by_operation(ctx, "differential_abundance").parameters["feature_count"] == 2
        @test _node_by_operation(ctx, "dia_like_quantification").parameters["feature_count"] == 2
        @test _node_by_operation(ctx, "protein_inference_top3").parameters["protein_count"] == 2
    end

    @testset "Enrichment" begin
        ctx = BioToolkit.ProvenanceContext()
        terms = BioToolkit.EnrichmentTerm[
            BioToolkit.EnrichmentTerm("GO:1", "term1", "GO", ["TP53", "BRCA1"], String[]),
            BioToolkit.EnrichmentTerm("GO:2", "term2", "GO", ["BRCA1"], ["GO:1"]),
        ]

        db = BioToolkit.build_annotation_database(terms; prov_ctx=ctx)
        results = BioToolkit.enrichment_test(["TP53", "BRCA1"], db; prov_ctx=ctx)
        ranked = [("TP53", 2.0), ("BRCA1", 1.0), ("MDM2", -1.0)]
        gsea = BioToolkit.gsea_preranked(ranked, db; prov_ctx=ctx, n_permutations=10, min_size=1, max_size=10)
        gsva = BioToolkit.gsva_score([1.0 2.0; 3.0 4.0], Dict("set1" => ["TP53", "BRCA1"]), ["TP53", "BRCA1"]; prov_ctx=ctx, min_size=1)
        propagated = BioToolkit.network_propagation(Dict("TP53" => 1.0), [0.0 1.0; 1.0 0.0], ["TP53", "BRCA1"]; prov_ctx=ctx)

        @test !isempty(results)
        @test !isempty(gsea)
        @test gsva isa BioToolkit.GSVAResult
        @test nrow(propagated) == 2

        @test _node_by_operation(ctx, "build_annotation_database").parameters["term_count"] == 2
        @test _node_by_operation(ctx, "enrichment_test").parameters["query_count"] == 2
        @test _node_by_operation(ctx, "fgsea_like").parameters["result_count"] >= 0
        @test _node_by_operation(ctx, "gsva_score").parameters["set_count"] == 1
        @test _node_by_operation(ctx, "network_propagation").parameters["gene_count"] == 2
    end

    @testset "Clinical" begin
        ctx = BioToolkit.ProvenanceContext()

        km = BioToolkit.kaplan_meier([1.0, 2.0, 3.0], [1, 0, 1]; prov_ctx=ctx)
        lr = BioToolkit.logrank_test([1.0, 2.0, 3.0, 4.0], [1, 1, 0, 1], ["A", "B", "A", "B"]; prov_ctx=ctx)
        pgx = BioToolkit.pharmacogenomics_star_alleles(DataFrame(gene=["CYP2D6"], variant=["rs1065852"]); prov_ctx=ctx)
        pheno = BioToolkit.cpic_metabolizer_phenotype(["*1/*4"]; prov_ctx=ctx)
        trial = BioToolkit.trial_suitability_scores(DataFrame(age=[55], ecog=[1]); prov_ctx=ctx)

        @test km isa BioToolkit.KaplanMeierResult
        @test lr.pvalue >= 0
        @test pgx.star_allele == ["*10"]
        @test pheno.phenotype == ["intermediate_metabolizer"]
        @test trial.eligible == [true]

        @test _node_by_operation(ctx, "kaplan_meier").parameters["n"] == 3
        @test _node_by_operation(ctx, "logrank_test").parameters["group_count"] == 2
        @test _node_by_operation(ctx, "pharmacogenomics_star_alleles").parameters["row_count"] == 1
        @test _node_by_operation(ctx, "cpic_metabolizer_phenotype").parameters["row_count"] == 1
        @test _node_by_operation(ctx, "trial_suitability_scores").parameters["row_count"] == 1
    end

    @testset "Explicit result provenance" begin
        distances = [0.0 1.0 2.0; 1.0 0.0 1.5; 2.0 1.5 0.0]
        pcoa = BioToolkit.pcoa(distances; dimensions=2)
        @test pcoa.provenance.source == "Microbiome/pcoa"
        @test pcoa.provenance.parameters.dimensions == 2

        nmds = BioToolkit.nmds(distances; dimensions=2, n_starts=2, maxiters=5, multi_thread=false)
        @test nmds.provenance.source == "Microbiome/nmds"
        @test nmds.provenance.parameters.n_starts == 2
        @test nmds.provenance.status in (:ok, :warn)

        assays = [[1.0 2.0; 3.0 4.0], [0.5 1.5; 2.5 3.5]]
        mofa = BioToolkit.multi_omics_factor_analysis(assays; n_factors=2)
        @test mofa.provenance.source == "SystemsBio/multi_omics_factor_analysis"
        @test mofa.provenance.parameters.assay_count == 2

        counts = [5 2 1; 3 4 2]
        sce = BioToolkit.SingleCellExperiment(counts, ["G1", "G2"], ["C1", "C2", "C3"])
        BioToolkit.attach_velocity_layers!(sce; spliced=counts, unspliced=counts .+ 1)
        velocity = BioToolkit.calculate_rna_velocity(sce; normalize=false)
        @test velocity.provenance.source == "SingleCell/calculate_rna_velocity"
        @test velocity.provenance.parameters.gene_count == 2

        dyn = BioToolkit.calculate_dynamical_rna_velocity(sce; normalize=false)
        @test dyn.provenance.source == "SingleCell/calculate_dynamical_rna_velocity"
        @test "latent time inferred from calculate_rna_velocity" in dyn.provenance.fallbacks

        wnn = BioToolkit.weighted_nearest_neighbors(Dict("rna" => sce, "adt" => sce); k=1, n_components=1)
        @test wnn.provenance.source == "SingleCell/weighted_nearest_neighbors"
        @test wnn.provenance.parameters.modality_count == 2

        perturb = BioToolkit.predict_perturbation(sce, "G1"; normalize=false, top_n=1)
        @test perturb.provenance.source == "SingleCell/predict_perturbation"
        @test perturb.provenance.parameters.target_gene == "G1"
    end

    @testset "Metadata provenance" begin
        coverage = Dict("chr1" => BioToolkit.SparseCoverageVector([1, 10, 20], [0, 6, 0], 20))
        peaks = BioToolkit.call_peaks(coverage; pvalue_threshold=1.0, min_depth=1)
        peak_prov = BioToolkit.metadata_provenance(peaks.metadata)
        @test peak_prov !== nothing
        @test peak_prov.source == "Epigenetics/call_peaks"

        calls = [
            BioToolkit.MethylationCall("chr1", 1, "s1", true),
            BioToolkit.MethylationCall("chr1", 2, "s1", false),
            BioToolkit.MethylationCall("chr1", 1, "s2", true),
        ]
        methylation = BioToolkit.bin_methylation(calls; bin_width=10)
        methylation_prov = BioToolkit.metadata_provenance(methylation.metadata)
        @test methylation_prov !== nothing
        @test methylation_prov.source == "Epigenetics/bin_methylation"

        alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite(BioToolkit.DNASeq("AAAA"); identifier="a"),
            BioToolkit.SeqRecordLite(BioToolkit.DNASeq("AAAT"); identifier="b"),
            BioToolkit.SeqRecordLite(BioToolkit.DNASeq("AATT"); identifier="c"),
        ])
        tree = BioToolkit.bootstrap_consensus_tree(alignment; replicates=2)
        tree_prov = BioToolkit.metadata_provenance(tree.metadata)
        @test tree_prov !== nothing
        @test tree_prov.source == "bootstrap_consensus_tree"
    end
end
