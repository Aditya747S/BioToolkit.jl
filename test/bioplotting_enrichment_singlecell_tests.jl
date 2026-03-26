@testset "BioPlotting" begin
    results = [
        BioToolkit.DEResult("gene1", 10.0, 2.5, 0.5, 12.0, 1e-6, 1e-4),
        BioToolkit.DEResult("gene2", 12.0, -2.1, 0.6, 11.0, 2e-4, 0.01),
        BioToolkit.DEResult("gene3", 8.0, 0.2, 0.4, 1.0, 0.8, 0.9),
    ]

    points, summary = BioToolkit.volcano_data(results; lfc_cutoff=1.0, fdr_cutoff=0.05, label_top=1)
    @test length(points) == 3
    @test summary.significant == 2
    @test any(point -> point.labeled, points)

    ma_points, ma_summary = BioToolkit.ma_data(results)
    @test length(ma_points) == 3
    @test ma_summary.significant == 2

    heatmap = BioToolkit.clustered_heatmap([1.0 2.0 3.0; 4.0 5.0 6.0]; row_labels=["r1", "r2"], column_labels=["c1", "c2", "c3"])
    @test size(heatmap.matrix) == (2, 3)
    @test length(heatmap.row_order) == 2
    @test length(heatmap.column_order) == 3
end

@testset "Enrichment" begin
    database = BioToolkit.build_annotation_database([
        BioToolkit.EnrichmentTerm("GO:0001", "Signal transduction", "GO", ["gene1", "gene2"], String[]),
        BioToolkit.EnrichmentTerm("GO:0002", "Cell cycle", "GO", ["gene2", "gene3"], ["GO:0001"]),
        BioToolkit.EnrichmentTerm("KEGG:0001", "Metabolism", "KEGG", ["gene3", "gene4"], String[]),
    ])

    builtin = BioToolkit.builtin_annotation_database()
    @test !isempty(builtin.terms)
    @test any(term -> term.namespace == "GO", values(builtin.terms))
    @test any(term -> term.namespace == "KEGG", values(builtin.terms))
    @test builtin.terms["GO:0006281"].parents == ["GO:0006974"]
    @test builtin.terms["GO:0006302"].parents == ["GO:0006281"]
    @test builtin.terms["GO:0000082"].parents == ["GO:0007049"]
    @test builtin.terms["KEGG:hsa05200"].parents == ["KEGG:hsa00001"]
    @test builtin.terms["KEGG:hsa04110"].parents == ["KEGG:hsa05200"]
    @test builtin.terms["GO:0006302"].genes ⊆ builtin.terms["GO:0006281"].genes
    @test builtin.terms["GO:0006281"].genes ⊆ builtin.terms["GO:0006974"].genes
    @test builtin.terms["GO:0007049"].genes ⊆ builtin.terms["GO:0008150"].genes
    @test builtin.terms["KEGG:hsa04110"].genes ⊆ builtin.terms["KEGG:hsa05200"].genes
    @test builtin.terms["KEGG:hsa05200"].genes ⊆ builtin.terms["KEGG:hsa00001"].genes

    builtin_go = BioToolkit.go_enrichment(["BRCA1", "BRCA2", "RAD51"], builtin)
    @test !isempty(builtin_go)
    @test any(result -> result.term_id == "GO:0006302", builtin_go)
    @test any(result -> result.term_id == "GO:0006281", builtin_go)

    mapped = BioToolkit.map_ids(database.mapper, ["gene1", "gene5"])
    @test mapped == ["gene1", "gene5"]

    results = BioToolkit.enrichment_test(["gene1", "gene2"], database)
    @test !isempty(results)
    @test results[1].term_id in ("GO:0001", "GO:0002")
    @test results[1].padj <= 1.0

    go_results = BioToolkit.go_enrichment(["gene1", "gene2"], database)
    @test all(result -> result.namespace == "GO", go_results)

    kegg_results = BioToolkit.kegg_enrichment(["gene3"], database)
    @test all(result -> result.namespace == "KEGG", kegg_results)

    dotplot_result = BioToolkit.dotplot(results; top_n=2)
    @test length(dotplot_result.results) <= 2
end

@testset "SingleCell" begin
    counts = [10 0 1 2; 12 1 0 3; 0 8 2 0; 1 1 7 9]
    experiment = BioToolkit.SingleCellExperiment(counts, ["gene1", "gene2", "gene3", "gene4"], ["cell1", "cell2", "cell3", "cell4"])

    matrix = BioToolkit.normalize_counts(experiment)
    @test size(matrix) == (4, 4)

    transformed, genes, cells = BioToolkit.sctransform(experiment; min_total=0)
    @test size(transformed, 2) == 4
    @test length(genes) == size(transformed, 1)
    @test length(cells) == 4

    pca = BioToolkit.run_pca(experiment; normalized=matrix, n_components=2)
    @test size(pca, 2) == 2
    @test haskey(experiment.reductions, "pca")

    umap = BioToolkit.run_umap(experiment; embedding=pca)
    @test size(umap, 2) <= 2

    labels = BioToolkit.cluster_cells(experiment; embedding=pca, method=:graph, k=2, random_seed=1)
    @test length(labels) == 4
    @test sum(values(BioToolkit.summarize_clusters(labels))) == 4

    marker_labels = BioToolkit.cluster_cells(experiment; embedding=pca, method=:kmeans, n_clusters=2, random_seed=1)
    @test length(unique(marker_labels)) == 2

    markers = BioToolkit.find_cluster_markers(experiment, marker_labels; cluster_id=marker_labels[1], min_total=0, shrink=false)
    @test !isempty(markers)
    @test all(result -> result isa BioToolkit.DEResult, markers)

    marker_summary = BioToolkit.cluster_marker_summary(experiment, marker_labels; min_total=0, shrink=false, top_n=2)
    @test !isempty(marker_summary)
    @test all(value -> isa(value, Vector{BioToolkit.DEResult}), values(marker_summary))
    @test all(length(value) <= 2 for value in values(marker_summary))
end